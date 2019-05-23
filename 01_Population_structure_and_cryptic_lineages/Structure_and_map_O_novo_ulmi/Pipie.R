pipie=function (x, long, lat, labels = names(x), edges = 50, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
    init.usr=par("usr")
    init.pin=par("pin")
x0=init.usr[1]
xmax=init.usr[2]
y0=init.usr[3]
ymax=init.usr[4]

#print(init.usr)
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("'x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(1L:length(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
#   plot.new()

    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)

      if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    plot.window(xlim, ylim, "", asp = 1)
new.usr=par("usr")
X0=new.usr[1]
Xmax=new.usr[2]
Y0=new.usr[3]
Ymax=new.usr[4]

long=((long-x0)*(Xmax-X0)/(xmax-x0))+X0
lat=((lat-y0)*(Ymax-Y0)/(ymax-y0))+Y0

#print(par("usr"))
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
                "lavender", "cornsilk")
        else par("fg")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }


    for (i in 1L:nx) {
#   for (i in 1L:1) {

        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
#	print(paste("P$x =", c(P$x,0)))
#	print(paste("P$y =", c(P$y,0)))
        polygon(c(P$x+long, long), c(P$y+lat, lat), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
        P <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
            text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
                adj = ifelse(P$x < 0, 1, 0), ...)
        }
    }
    title(main = main, ...)
    invisible(NULL)
    par(usr=init.usr)
    par(pin=init.pin)

}

