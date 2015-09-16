smilyface <-
function (center, r, mood = c("happy", "sad")) {
    mood <- match.arg(mood)
    getYmult <- function () {
        if (dev.cur() == 1) {
            warning("No graphics device open.")
            ymult <- 1
        }
        else {
            xyasp <- par("pin")
            xycr <- diff(par("usr"))[c(1, 3)]
            ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
        }
        return(ymult)
    }
    draw.circle <- function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, 
                             lwd = 1) {
        xylim <- par("usr")
        plotdim <- par("pin")
        ymult <- getYmult()
        angle.inc <- 2 * pi/nv
        angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
        if (length(col) < length(radius)) 
            col <- rep(col, length.out = length(radius))
        for (circle in 1:length(radius)) {
            xv <- cos(angles) * radius[circle] + x
            yv <- sin(angles) * radius[circle] * ymult + y
            polygon(xv, yv, border = border, col = col[circle], lty = lty, 
                    lwd = lwd)
        }
        invisible(list(x = xv, y = yv))
    }
    draw.arc <- function (x = 1, y = NULL, radius = 1, angle1 = deg1 * pi/180, 
                          angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 0.05, col = NA, 
                          lwd = NA, ...) {
        if (all(is.na(col))) 
            col <- par("col")
        if (all(is.na(lwd))) 
            lwd <- par("lwd")
        xylim <- par("usr")
        ymult <- getYmult()
        devunits <- dev.size("px")
        draw.arc.0 <- function(x, y, radius, angle1, angle2, n, col, 
                               lwd, ...) {
            delta.angle <- (angle2 - angle1)
            if (n != as.integer(n)) 
                n <- as.integer(1 + delta.angle/n)
            delta.angle <- delta.angle/n
            angleS <- angle1 + seq(0, length = n) * delta.angle
            angleE <- c(angleS[-1], angle2)
            if (n > 1) {
                half.lwd.user <- (lwd/2) * (xylim[2] - xylim[1])/devunits[1]
                adj.angle = delta.angle * half.lwd.user/(2 * (radius + 
                                                                  half.lwd.user))
                angleS[2:n] = angleS[2:n] - adj.angle
                angleE[1:(n - 1)] = angleE[1:(n - 1)] + adj.angle
            }
            p1x <- x + radius * cos(angleS)
            p1y <- y + radius * sin(angleS) * ymult
            p2x <- x + radius * cos(angleE)
            p2y <- y + radius * sin(angleE) * ymult
            segments(p1x, p1y, p2x, p2y, col = col, lwd = lwd, ...)
        }
        xy <- xy.coords(x, y)
        x <- xy$x
        y <- xy$y
        a1 <- pmin(angle1, angle2)
        a2 <- pmax(angle1, angle2)
        angle1 <- a1
        angle2 <- a2
        args <- data.frame(x, y, radius, angle1, angle2, n, col, 
                           lwd, stringsAsFactors = FALSE)
        for (i in 1:nrow(args)) do.call("draw.arc.0", c(args[i, ], 
                                                        ...))
        invisible(args)
    }
    if (mood == "happy") {
        draw.circle(center[1], center[2], radius = r, col = "green")
        points(center[1] - r/3, center[2] + r/3, pch = 16, cex = 0.5)
        points(center[1] + r/3, center[2] + r/3, pch = 16, cex = 0.5)
        draw.arc(center[1], center[2], radius = r/2, deg1 = 215, deg = 325)
    } else {
        draw.circle(center[1], center[2], radius = r, col = "red")
        points(center[1] - r/3, center[2] + r/3, pch = 16, cex = 0.5)
        points(center[1] + r/3, center[2] + r/3, pch = 16, cex = 0.5)
        draw.arc(center[1], center[2] - 4.5*r, radius = r/2, deg1 = 145, deg = 35)
    }
}
