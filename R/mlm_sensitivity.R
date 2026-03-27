#' Robustness diagnostics for cross-level interaction effects
#'
#' Assesses the stability of a cross-level interaction effect using two
#' MLM-appropriate diagnostics:
#'
#' 1. **ICC-shift robustness** -- how do the interaction SE and Johnson-Neyman
#'    boundary change if the intraclass correlation were different from what
#'    was observed? This is relevant because the ICC determines the effective
#'    sample size at level 2, which directly governs precision of cross-level
#'    interaction estimates.
#'
#' 2. **Leave-one-cluster-out (LOCO) stability** -- refit the model dropping
#'    one cluster at a time and track how the interaction coefficient moves.
#'    This is nonparametric, makes no distributional assumptions, and directly
#'    answers: "Is this finding driven by a small number of influential
#'    clusters?"
#'
#' @section Scope:
#' These are **robustness diagnostics**, not a full causal sensitivity
#' analysis. They do not quantify the strength of unmeasured confounding
#' needed to explain away the interaction -- that requires a level-2-aware
#' omitted variable bound that is currently under development as a separate
#' methodological contribution. See `vignette("robustness-diagnostics")` for
#' interpretation guidance.
#'
#' @param model An `lmerMod` object with a two-way interaction between `pred`
#'   and `modx`.
#' @param pred Character scalar. Focal predictor name.
#' @param modx Character scalar. Moderator name.
#' @param alpha Significance level. Default `0.05`.
#' @param icc_range Numeric vector of length 2. Range of ICC values to
#'   evaluate. Default `c(0.01, 0.40)`.
#' @param icc_grid Integer. Number of ICC values in the grid. Default `50`.
#' @param loco Logical. Run leave-one-cluster-out analysis? Default `TRUE`.
#'   Set to `FALSE` for large datasets where refitting is slow.
#' @param conf.level Confidence level. Default `0.95`.
#' @param verbose Logical. Print progress during LOCO refitting? Default
#'   `FALSE`.
#'
#' @return An object of class `mlm_sensitivity` with components:
#'   * `icc_shift`: data frame of interaction SE, t, p, significance, and
#'     approximate JN boundary across the ICC grid.
#'   * `loco`: data frame with one row per cluster giving the interaction
#'     coefficient, SE, t, p, and Cook's-distance-style influence measure
#'     when that cluster is omitted. `NULL` if `loco = FALSE`.
#'   * `robustness_index`: proportion of the ICC range where the interaction
#'     remains significant.
#'   * `observed`: list of observed model statistics.
#'   * Metadata: `pred`, `modx`, `alpha`, `icc_range`, `int_term`.
#'
#' @examples
#' # Use a small dataset for fast execution
#' set.seed(42)
#' n_j <- 20; n_i <- 10
#' dat_small <- data.frame(
#'   y    = rnorm(n_j * n_i),
#'   x    = rnorm(n_j * n_i),
#'   m    = rep(rnorm(n_j), each = n_i),
#'   grp  = factor(rep(seq_len(n_j), each = n_i))
#' )
#' dat_small$y <- dat_small$y + 0.5 * dat_small$x * dat_small$m
#' mod_small <- lme4::lmer(y ~ x * m + (1 + x | grp), data = dat_small,
#'                         control = lme4::lmerControl(optimizer = "bobyqa"))
#'
#' # ICC-shift only (fast)
#' sens <- mlm_sensitivity(mod_small, pred = "x", modx = "m", loco = FALSE)
#' print(sens)
#'
#' # Full diagnostics including LOCO (20 clusters - fast)
#' \donttest{
#' sens_full <- mlm_sensitivity(mod_small, pred = "x", modx = "m")
#' plot(sens_full)
#' }
#'
#' @export
mlm_sensitivity <- function(model,
                             pred,
                             modx,
                             alpha      = 0.05,
                             icc_range  = c(0.01, 0.40),
                             icc_grid   = 50L,
                             loco       = TRUE,
                             conf.level = 0.95,
                             verbose    = FALSE) {

  .check_lmer(model)
  .validate_terms(model, pred, modx)

  int_term <- .get_interaction_term(model, pred, modx)
  fe       <- lme4::fixef(model)
  vcv      <- .extract_vcov(model)
  df_r     <- get_residual_df(model)

  b_int  <- fe[int_term]
  se_int <- sqrt(vcv[int_term, int_term])
  t_int  <- b_int / se_int
  p_int  <- 2 * stats::pt(abs(t_int), df = df_r, lower.tail = FALSE)

  # Observed JN bounds
  jn_out    <- mlm_jn(model, pred = pred, modx = modx, alpha = alpha)
  jn_bounds <- jn_out$jn_bounds

  # Observed ICC
  vc         <- lme4::VarCorr(model)
  grp_name   <- names(lme4::getME(model, "flist"))[1]
  tau00      <- as.numeric(vc[[grp_name]][1, 1])
  sigma2     <- as.numeric(attr(vc, "sc")^2)
  icc_obs    <- tau00 / (tau00 + sigma2)

  mf         <- model@frame
  n_total    <- nrow(mf)
  n_clusters <- length(unique(mf[[grp_name]]))
  n_avg      <- n_total / n_clusters

  # ------ ICC-shift analysis ------------------------------------------------------------------------------------------------------------------------------------------------------------
  icc_seq  <- seq(icc_range[1], icc_range[2], length.out = icc_grid)
  deff_obs <- 1 + (n_avg - 1) * icc_obs

  icc_rows <- lapply(icc_seq, function(rho) {
    deff_adj <- 1 + (n_avg - 1) * rho
    se_adj   <- se_int * sqrt(deff_adj / deff_obs)
    t_adj    <- b_int / se_adj
    p_adj    <- 2 * stats::pt(abs(t_adj), df = df_r, lower.tail = FALSE)
    sig      <- p_adj < alpha
    jn_adj   <- if (!all(is.na(jn_bounds))) jn_bounds[1] * (se_adj / se_int) else NA_real_

    data.frame(
      icc          = rho,
      deff         = deff_adj,
      se_int_adj   = se_adj,
      t_int_adj    = t_adj,
      p_int_adj    = p_adj,
      sig          = sig,
      jn_bound_adj = jn_adj,
      stringsAsFactors = FALSE
    )
  })
  icc_df <- do.call(rbind, icc_rows)
  rownames(icc_df) <- NULL

  robustness_index <- mean(icc_df$sig, na.rm = TRUE)

  # ------ Leave-one-cluster-out (LOCO) ------------------------------------------------------------------------------------------------------------------------------
  loco_df <- NULL

  if (loco) {
    cluster_ids <- unique(mf[[grp_name]])
    formula_mod <- stats::formula(model)

    if (verbose) cat("Running LOCO across", length(cluster_ids), "clusters...\n")

    loco_rows <- lapply(seq_along(cluster_ids), function(i) {
      cid <- cluster_ids[i]
      if (verbose && i %% 10 == 0) cat("  Cluster", i, "of", length(cluster_ids), "\n")

      mf_sub <- mf[mf[[grp_name]] != cid, , drop = FALSE]

      fit_sub <- tryCatch(
        suppressMessages(
          lme4::lmer(formula_mod, data = mf_sub, REML = FALSE,
                     control = lme4::lmerControl(optimizer = "bobyqa"))
        ),
        error = function(e) NULL
      )

      if (is.null(fit_sub)) {
        return(data.frame(
          cluster     = as.character(cid),
          b_int       = NA_real_,
          se_int      = NA_real_,
          t_int       = NA_real_,
          p_int       = NA_real_,
          b_change    = NA_real_,
          pct_change  = NA_real_,
          sig         = NA,
          stringsAsFactors = FALSE
        ))
      }

      fe_sub   <- lme4::fixef(fit_sub)
      vcv_sub  <- as.matrix(stats::vcov(fit_sub))

      # Handle interaction term name (order may vary)
      int_sub <- if (int_term %in% names(fe_sub)) {
        int_term
      } else {
        paste0(modx, ":", pred)
      }

      b_sub  <- fe_sub[int_sub]
      se_sub <- sqrt(vcv_sub[int_sub, int_sub])
      df_sub <- nrow(mf_sub) - length(fe_sub)
      t_sub  <- b_sub / se_sub
      p_sub  <- 2 * stats::pt(abs(t_sub), df = df_sub, lower.tail = FALSE)

      data.frame(
        cluster    = as.character(cid),
        b_int      = b_sub,
        se_int     = se_sub,
        t_int      = t_sub,
        p_int      = p_sub,
        b_change   = b_sub - b_int,
        pct_change = 100 * (b_sub - b_int) / abs(b_int),
        sig        = p_sub < alpha,
        stringsAsFactors = FALSE
      )
    })

    loco_df <- do.call(rbind, loco_rows)
    rownames(loco_df) <- NULL

    # Cook's-distance-style influence: standardised change in b_int
    sd_change      <- stats::sd(loco_df$b_change, na.rm = TRUE)
    loco_df$influence <- abs(loco_df$b_change) / (sd_change + .Machine$double.eps)
  }

  structure(
    list(
      observed = list(
        b_int      = b_int,
        se_int     = se_int,
        t_int      = t_int,
        p_int      = p_int,
        jn_bounds  = jn_bounds,
        icc_obs    = icc_obs,
        n_total    = n_total,
        n_clusters = n_clusters,
        n_avg      = n_avg
      ),
      icc_shift        = icc_df,
      loco             = loco_df,
      robustness_index = robustness_index,
      pred             = pred,
      modx             = modx,
      alpha            = alpha,
      icc_range        = icc_range,
      conf.level       = conf.level,
      grp_name         = grp_name,
      int_term         = int_term
    ),
    class = "mlm_sensitivity"
  )
}

#' @export
print.mlm_sensitivity <- function(x, digits = 3, ...) {

  cat("\n========================================\n")
  cat("  Robustness Diagnostics\n")
  cat("  (Cross-Level Interaction)\n")
  cat("========================================\n")
  cat("Focal predictor :", x$pred, "\n")
  cat("Moderator       :", x$modx, "\n")
  cat("Interaction term:", x$int_term, "\n\n")

  obs <- x$observed
  cat("--- Observed interaction ---\n")
  cat(sprintf("  b      = %7.3f\n", obs$b_int))
  cat(sprintf("  SE     = %7.3f\n", obs$se_int))
  cat(sprintf("  t      = %7.3f\n", obs$t_int))
  cat(sprintf("  p      = %s\n",    format_pval(obs$p_int)))
  cat(sprintf("  ICC    = %7.3f\n", obs$icc_obs))
  cat(sprintf("  J      = %d clusters\n", obs$n_clusters))
  cat(sprintf("  n_avg  = %.1f per cluster\n", obs$n_avg))

  if (!all(is.na(obs$jn_bounds))) {
    cat(sprintf("  JN boundary: %s\n",
                paste(round(obs$jn_bounds, 3), collapse = ", ")))
  }

  cat("\n--- ICC-shift robustness ---\n")
  cat(sprintf("  ICC range: [%.2f, %.2f]\n",
              x$icc_range[1], x$icc_range[2]))
  cat(sprintf("  Robustness index: %.1f%%\n", x$robustness_index * 100))
  cat(sprintf("  (interaction significant across %.1f%% of ICC range)\n",
              x$robustness_index * 100))
  se_range <- range(x$icc_shift$se_int_adj, na.rm = TRUE)
  cat(sprintf("  SE ranges from %.3f to %.3f across ICC grid\n",
              se_range[1], se_range[2]))

  if (!is.null(x$loco)) {
    loco <- x$loco
    n_sig     <- sum(loco$sig, na.rm = TRUE)
    n_valid   <- sum(!is.na(loco$sig))
    pct_sig   <- 100 * n_sig / n_valid
    b_range   <- range(loco$b_int, na.rm = TRUE)
    top_inf   <- loco[order(-loco$influence), ][1:min(3, nrow(loco)), ]

    cat("\n--- Leave-one-cluster-out (LOCO) ---\n")
    cat(sprintf("  Interaction significant in %d / %d fits (%.1f%%)\n",
                n_sig, n_valid, pct_sig))
    cat(sprintf("  b ranges from %.3f to %.3f\n", b_range[1], b_range[2]))
    cat(sprintf("  Max |change| = %.3f (%.1f%% of original b)\n",
                max(abs(loco$b_change), na.rm = TRUE),
                max(abs(loco$pct_change), na.rm = TRUE)))

    cat("\n  Most influential clusters:\n")
    cat(sprintf("  %-12s  %8s  %10s  %10s\n",
                "Cluster", "b_int", "b_change", "% change"))
    for (i in seq_len(nrow(top_inf))) {
      cat(sprintf("  %-12s  %8.3f  %10.3f  %10.1f\n",
                  top_inf$cluster[i], top_inf$b_int[i],
                  top_inf$b_change[i], top_inf$pct_change[i]))
    }
  }

  if (x$robustness_index > 0.90) {
    cat("\nOVERALL: Interaction is stable across the ICC range tested.\n")
  } else if (x$robustness_index > 0.60) {
    cat("\nOVERALL: Interaction is moderately sensitive to ICC assumptions.\n")
  } else {
    cat("\nOVERALL: Interaction is fragile -- sensitive to ICC assumptions.\n")
  }

  cat("\nNOTE: These are robustness diagnostics, not a full causal\n")
  cat("sensitivity analysis. Level-2 omitted variable bounds for\n")
  cat("cross-level interactions are under methodological development.\n\n")

  invisible(x)
}

#' Plot robustness diagnostics for a cross-level interaction
#'
#' Produces up to three panels: (1) interaction SE across the ICC range,
#' (2) JN boundary shift across the ICC range, and (3) LOCO coefficient
#' stability plot showing the interaction estimate when each cluster is
#' omitted, with influential clusters flagged.
#'
#' @param x An `mlm_sensitivity` object.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot.mlm_sensitivity <- function(x, ...) {

  df      <- x$icc_shift
  obs_icc <- x$observed$icc_obs
  alpha   <- x$alpha
  b_obs   <- x$observed$b_int

  # ------ Panel 1: ICC shift - SE ---------------------------------------------------------------------------------------------------------------------------------------------
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = icc, y = se_int_adj)) +
    ggplot2::geom_line(ggplot2::aes(color = sig), linewidth = 1.1) +
    ggplot2::geom_vline(xintercept = obs_icc, linetype = "dashed",
                        color = "grey30", linewidth = 0.8) +
    ggplot2::annotate("text", x = obs_icc, y = max(df$se_int_adj, na.rm = TRUE),
                      label = paste0("Observed\nICC = ", round(obs_icc, 3)),
                      hjust = -0.1, size = 3, color = "grey30") +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
      labels = c("TRUE" = paste0("p < ", alpha),
                 "FALSE" = paste0("p >= ", alpha)),
      name = NULL
    ) +
    ggplot2::labs(x = "ICC", y = "SE of interaction",
                  title = "ICC Sensitivity: Interaction SE") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position  = "bottom",
                   panel.grid.major = ggplot2::element_line(color = "grey92"))

  # ------ Panel 2: ICC shift - JN boundary ---------------------------------------------------------------------------------------------------------------
  has_jn <- !all(is.na(df$jn_bound_adj))

  p2 <- if (has_jn) {
    ggplot2::ggplot(df, ggplot2::aes(x = icc, y = jn_bound_adj)) +
      ggplot2::geom_line(color = "#7B2D8B", linewidth = 1.1) +
      ggplot2::geom_vline(xintercept = obs_icc, linetype = "dashed",
                          color = "grey30", linewidth = 0.8) +
      ggplot2::labs(x = "ICC",
                    y = paste0("JN boundary (", x$modx, ")"),
                    title = "ICC Sensitivity: JN Boundary") +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = "grey92"))
  } else NULL

  # ------ Panel 3: LOCO stability ---------------------------------------------------------------------------------------------------------------------------------------------
  p3 <- NULL

  if (!is.null(x$loco)) {
    loco <- x$loco
    loco <- loco[!is.na(loco$b_int), ]

    # Flag top 3 most influential
    loco$label <- ""
    top3 <- order(-loco$influence)[1:min(3, nrow(loco))]
    loco$label[top3] <- loco$cluster[top3]
    loco$influential <- loco$influence > stats::quantile(loco$influence,
                                                          0.90, na.rm = TRUE)

    p3 <- ggplot2::ggplot(loco,
             ggplot2::aes(x = stats::reorder(cluster, b_int), y = b_int)) +
      # Reference line: full-model estimate
      ggplot2::geom_hline(yintercept = b_obs, linetype = "dashed",
                          color = "#2166AC", linewidth = 0.8) +
      ggplot2::annotate("text", x = 1, y = b_obs,
                        label = paste0("Full model\nb = ", round(b_obs, 3)),
                        hjust = 0, vjust = -0.4, size = 3, color = "#2166AC") +
      # LOCO estimates
      ggplot2::geom_point(ggplot2::aes(color = influential),
                          size = 1.5, alpha = 0.7) +
      # Label top 3
      ggplot2::geom_text(ggplot2::aes(label = label),
                         hjust = -0.2, size = 2.8, color = "#D6604D") +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "grey50", "TRUE" = "#D6604D"),
        labels = c("FALSE" = "Typical", "TRUE" = "Influential (top 10%)"),
        name   = NULL
      ) +
      ggplot2::scale_x_discrete(labels = NULL, breaks = NULL) +
      ggplot2::labs(
        x        = paste0("Clusters (n = ", nrow(loco), "), ordered by b"),
        y        = paste0("b (", x$int_term, ")"),
        title    = "LOCO: Interaction Stability Across Clusters",
        subtitle = paste0(
          sum(loco$sig, na.rm = TRUE), "/", nrow(loco),
          " fits remain significant  |  ",
          "b range: [", round(min(loco$b_int, na.rm = TRUE), 3),
          ", ", round(max(loco$b_int, na.rm = TRUE), 3), "]"
        )
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position  = "bottom",
        panel.grid.major = ggplot2::element_line(color = "grey92"),
        plot.subtitle    = ggplot2::element_text(color = "grey40", size = 9)
      )
  }

  # ------ Combine ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("Install 'patchwork' for the full multi-panel layout.")
    return(if (!is.null(p3)) p3 else p1)
  }

  top_row <- if (!is.null(p2)) (p1 + p2) else p1

  result <- if (!is.null(p3)) {
    top_row / p3 +
      patchwork::plot_annotation(
        title    = "Robustness Diagnostics: Cross-Level Interaction",
        subtitle = paste0(x$int_term, "  |  b = ", round(b_obs, 3),
                          "  |  ICC = ", round(obs_icc, 3),
                          "  |  J = ", x$observed$n_clusters, " clusters")
      )
  } else {
    top_row
  }

  result
}
