lnl_tp_olsen <- function(param, X, y, wt, sum = TRUE, gradient = FALSE, hessian = FALSE,
                         left = 0, right = Inf, sample = "censored"){

    CENS <- sample == "censored"
    TRUNC <- ! CENS
    LEFT <- ! is.infinite(left) & is.infinite(right)
    RIGHT <- is.infinite(left) & ! is.infinite(right)
    TWO <- ! is.infinite(left) & ! is.infinite(right)
    if (TWO & TRUNC) stop("the two-limits truncated model is not supported")
    Ia <- y <= left
    Ib <- y >= right
    Io <- (y > left & y < right)

    K <- length(param) - 1
    N <- length(y)
    beta <- param[1:K]
    sig <- param[K + 1]

    bX <- as.numeric(X %*% beta)
    za <- (sig * left - bX)
    zb <- (sig * right - bX)
    ze <- (sig * y - bX)

    lnl_com <-  Io * ( log(sig) - 0.5 * log(2 * pi) - 1 / 2 * ze ^ 2)

    if (CENS){
        lnl_spec <- 0
        if (LEFT + TWO)  lnl_spec <-            Ia * pnorm(za, log.p = TRUE)
        if (RIGHT + TWO) lnl_spec <- lnl_spec + Ib * pnorm(- zb, log.p = TRUE)
    }
    if (TRUNC){
        if (LEFT)  lnl_spec <- - Io * pnorm(- za, log.p = TRUE)
        if (RIGHT) lnl_spec <- - Io * pnorm(  zb, log.p = TRUE)
        if (TWO)   lnl_spec <- - Io * log(pnorm(zb) - pnorm(za))
    }
    lnl <- lnl_com + lnl_spec
    if (sum) lnl <- sum(wt * lnl)
    
    if (gradient){
        g_com_beta <- Io * ze
        g_com_sig <-  Io * (1 / sig  - ze * y)
        if (CENS){
            g_spec_beta <- g_spec_sig <- 0
            if (LEFT | TWO){
                g_spec_beta <-  - Ia * mills(za)
                g_spec_sig <-     Ia * mills(za) * left
            }
            if(RIGHT | TWO){
                g_spec_beta <-   g_spec_beta +  Ib * mills(- zb)
                g_spec_sig <-    g_spec_sig  -  Ib * mills(-zb) * right
            }
        }
        if (TRUNC){
            if (LEFT){
                g_spec_beta <- -  Io * mills(- za)
                g_spec_sig <-     Io * mills(- za) * left
            }
            if (RIGHT){
                g_spec_beta <-    Io * mills(zb)
                g_spec_sig <-   - Io * mills(zb) * right
            }
            if (TWO){
                Dphi <- pnorm(zb) - pnorm(za)
                g_spec_beta <-   Io * (mills(zb) * pnorm(zb) - mills(za) * pnorm(za)) / Dphi
                g_spec_sig <-  - Io * (mills(zb) * pnorm(zb) * right - mills(za) * pnorm(za) * left) / Dphi
            }
        }
        grad <- wt * cbind((g_com_beta +  g_spec_beta) * X,
                           g_com_sig + g_spec_sig)
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    
    if (hessian){
        h_com_beta_beta <- - Io
        h_com_beta_sig <-  Io * y
        h_com_sig_sig <- - Io * (1 / sig ^ 2 + y ^ 2)
        if (CENS){
            h_spec_beta_beta <- h_spec_beta_sig <- h_spec_sig_sig <- 0
            if (LEFT | TWO){
                h_spec_beta_beta <-  Ia * dmills(za)
                h_spec_beta_sig <-   - Ia * dmills(za) * left
                h_spec_sig_sig <-  Ia * dmills(za) * left ^ 2
            }
            if (RIGHT | TWO){
                h_spec_beta_beta <- h_spec_beta_beta + Ib * dmills(-zb)
                h_spec_beta_sig <-  h_spec_beta_sig  - Ib * dmills(- zb) * right
                h_spec_sig_sig <-   h_spec_sig_sig   + Ib * dmills(- zb) * right ^ 2
            }
        }
        if (TRUNC){
            if (LEFT){
                h_spec_beta_beta <- - Io * dmills( - za) 
                h_spec_beta_sig <-    Io * left * dmills(- za)
                h_spec_sig_sig <-   - Io * dmills(- za) * left ^ 2
            }
            if (RIGHT){
                h_spec_beta_beta <- - Io * dmills(zb)
                h_spec_beta_sig <-    Io * dmills(zb) * right
                h_spec_sig_sig <-   - Io * dmills(zb) * right ^ 2
            }
            if (TWO){
                DELTA <- pnorm(zb) - pnorm(za)
                A <- mills(zb) * pnorm(zb) - mills(za) * pnorm(za)
                B <- mills(zb) * pnorm(zb) * right - mills(za) * pnorm(za) * left
                MILLS <- A / DELTA
                A_sig <- - pnorm(zb) * (mills(zb) ^ 2 + dmills(zb)) * zb + pnorm(za) * (mills(za) ^ 2 + dmills(za)) * za
                B_sig <- - pnorm(zb) * zb * (dmills(zb) * zb + mills(zb) ^ 2 * zb + mills(zb)) +
                     pnorm(za) * za * (dmills(za) * za + mills(za) ^ 2 * za + mills(za))
                D_sig <- - mills(zb) * pnorm(zb) * zb + mills(za) * pnorm(za) * za
                g_spec_beta <- - Io * MILLS
                g_spec_sig <- - Io * (mills(zb) * pnorm(zb) * zb - mills(za) * pnorm(za) * za) / Dphi
                DD <-  (pnorm(zb) * (dmills(zb) + mills(zb) ^ 2) - pnorm(za) * (dmills(za) + mills(za) ^ 2)) / Dphi -
                    (mills(zb) * pnorm(zb) - mills(za) * pnorm(za)) ^ 2 / (pnorm(zb) - pnorm(za)) ^ 2
                h_spec_beta_beta <- Io * DD
                h_spec_beta_sig <-   Io * (MILLS - (A_sig - MILLS * D_sig) / DELTA)
                h_spec_sig_sig <- Io * (B / DELTA - (B_sig * DELTA - B * D_sig) / DELTA ^ 2)
            }
        }
        H_bb <- crossprod(wt * (h_com_beta_beta + h_spec_beta_beta) * X, X)
        H_bs <- apply(wt * (h_com_beta_sig + h_spec_beta_sig) * X, 2, sum)
        H_ss <- sum(wt * (h_com_sig_sig + h_spec_sig_sig))
        attr(lnl, "hessian") <- rbind(cbind(H_bb, H_bs),
                                      c(H_bs, H_ss))
      }
    lnl
}
    


