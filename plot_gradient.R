library(ncdf4)


#d0 <- nc_open("/home/cgaillard/NOF/pop_1_0.nc")
#d1 <- nc_open("/home/cgaillard/WF/pop_1_1.nc")
#d0 <- nc_open("/home/cgaillard/Test-reference-NOFire/pop_1_0.nc")
d0 <- nc_open("pop_1_0.nc")

#d1 <- nc_open("/home/cgaillard/Test-decomposition-rate/pop_1_1.nc")
#d1 <- nc_open("/home/cgaillard/Test-no-fire-100years/pop_1_1.nc")
#d1 <- nc_open("/home/cgaillard/Test-intensity-300/pop_1_1.nc")
#d1 <- nc_open("/home/cgaillard/Test-reference/pop_1_1.nc")
#d1 <- nc_open("/home/cgaillard/Test-return-interval-7years/pop_1_1.nc")
d1 <- nc_open("pop_1_1.nc")

#mv transect_lat.pdf transect_lat_return-interval-7years.pdf
#mv transect_lat.pdf transect_lat_Test-reference.pdf
#mv transect_lat.pdf transect_lat_intensity-300.pdf
#mv transect_lat.pdf transect_lat_no-fire-100years.pdf
#mv transect_lat.pdf transect_lat_decomposition-rate.pdf


x_coo <- ncvar_get( d0, "lon" )
y_coo <- ncvar_get( d0, "lat" )
tt    <- ncvar_get( d0, "time" )
ll    <- 5*12

#sn0 <- ncvar_get( d0, "MeanStem", start=c( 1,1,2, length(tt) ), count=c( length(x_coo), length(y_coo), 1, 1 ) ) 
#sn1 <- ncvar_get( d1, "MeanStem", start=c( 1,1,2, length(tt) ), count=c( length(x_coo), length(y_coo), 1, 1 ) ) 
#bw0 <- ncvar_get( d0, "SumBWood", start=c( 1,1,2, length(tt) ), count=c( length(x_coo), length(y_coo), 1, 1 ) ) 
#bw1 <- ncvar_get( d1, "SumBWood", start=c( 1,1,2, length(tt) ), count=c( length(x_coo), length(y_coo), 1, 1 ) ) 
#bl0 <- ncvar_get( d0, "SumBLeaf", start=c( 1,1,3, length(tt) ), count=c( length(x_coo), length(y_coo), 1, 1 ) ) 
#bl1 <- ncvar_get( d1, "SumBLeaf", start=c( 1,1,3, length(tt) ), count=c( length(x_coo), length(y_coo), 1, 1 ) ) 
#vt0 <- ncvar_get( d0, "MeanVegType", start=c( 1,1, length(tt) ), count=c( length(x_coo), length(y_coo), 1 ) ) 
#vt1 <- ncvar_get( d1, "MeanVegType", start=c( 1,1, length(tt) ), count=c( length(x_coo), length(y_coo), 1 ) ) 

sn0t <- ( ncvar_get( d0, "MeanStem",    start=c( 1,1,2, length(tt)-ll ), count=c( length(x_coo), length(y_coo), 1, ll ) ) )
sn1t <- ( ncvar_get( d1, "MeanStem",    start=c( 1,1,2, length(tt)-ll ), count=c( length(x_coo), length(y_coo), 1, ll ) ) )
bw0t <- ( ncvar_get( d0, "SumBWood",    start=c( 1,1,2, length(tt)-ll ), count=c( length(x_coo), length(y_coo), 1, ll ) ) )
bw1t <- ( ncvar_get( d1, "SumBWood",    start=c( 1,1,2, length(tt)-ll ), count=c( length(x_coo), length(y_coo), 1, ll ) ) )
bl0t <- ( ncvar_get( d0, "SumBLeaf",    start=c( 1,1,3, length(tt)-ll ), count=c( length(x_coo), length(y_coo), 1, ll ) ) )
bl1t <- ( ncvar_get( d1, "SumBLeaf",    start=c( 1,1,3, length(tt)-ll ), count=c( length(x_coo), length(y_coo), 1, ll ) ) )
vt0t <- ( ncvar_get( d0, "MeanVegType", start=c( 1,1,   length(tt)-ll ), count=c( length(x_coo), length(y_coo),    ll ) ) )
vt1t <- ( ncvar_get( d1, "MeanVegType", start=c( 1,1,   length(tt)-ll ), count=c( length(x_coo), length(y_coo),    ll ) ) )
sn0 <- apply( sn0t, c(1,2), mean )
sn1 <- apply( sn1t, c(1,2), mean )
bw0 <- apply( bw0t, c(1,2), mean )
bw1 <- apply( bw1t, c(1,2), mean )
bl0 <- apply( bl0t, c(1,2), mean )
bl1 <- apply( bl1t, c(1,2), mean )
vt0 <- apply( vt0t, c(1,2), mean )
vt1 <- apply( vt1t, c(1,2), mean )


graphics.off()
pdf( "transect_lat.pdf", height=8, width=14 )
par( mfcol=c(2,2), mar=c(4,4,1,1) )

plot( -100, -100, xlim=range(y_coo), ylim=c(0,10), xlab="y coordinate", ylab="Mean stem number" )
for ( i in 1:2 ) points( y_coo, sn0[i,], col="steelblue", pch=15 )
for ( i in 1:2 ) points( y_coo, sn1[i,], col="darkgoldenrod", pch=15 )
lines( y_coo, apply( sn0, 2, mean, na.rm=T), lwd=2, col="steelblue" )
lines( y_coo, apply( sn1, 2, mean, na.rm=T), lwd=2, col="darkgoldenrod" )

#--------------------------------

plot( -100, -100, xlim=range(y_coo), ylim=c(0,1), xlab="y coordinate", ylab="Tree:grass ratio" )
for ( i in 1:2 ) points( y_coo, vt0[i,], col="steelblue", pch=15 )
for ( i in 1:2 ) points( y_coo, vt1[i,], col="darkgoldenrod", pch=15 )
lines( y_coo, apply( vt0, 2, mean, na.rm=T), lwd=2, col="steelblue" )
lines( y_coo, apply( vt1, 2, mean, na.rm=T), lwd=2, col="darkgoldenrod" )

#--------------------------------

plot( -100, -100, xlim=range(y_coo), ylim=c(0,400), xlab="y coordinate", ylab="Tree woody biomass" )
for ( i in 1:2 ) points( y_coo, bw0[i,], col="steelblue", pch=15,  )
for ( i in 1:2 ) points( y_coo, bw1[i,], col="darkgoldenrod", pch=15 )
lines( y_coo, apply( bw0, 2, mean, na.rm=T), lwd=2, col="steelblue" )
lines( y_coo, apply( bw1, 2, mean, na.rm=T), lwd=2, col="darkgoldenrod" )

#--------------------------------

plot( -100, -100, xlim=range(y_coo), ylim=c(0,2.5), xlab="y coordinate", ylab="Grass biomass" )
for ( i in 1:2 ) points( y_coo, bl0[i,], col="steelblue", pch=15 )
for ( i in 1:2 ) points( y_coo, bl1[i,], col="darkgoldenrod", pch=15 )
lines( y_coo, apply( bl0, 2, mean, na.rm=T), lwd=2, col="steelblue" )
lines( y_coo, apply( bl1, 2, mean, na.rm=T), lwd=2, col="darkgoldenrod" )

legend("topleft", c("No fire", "With fire"), col=c("steelblue", "darkgoldenrod"), lwd=2, bty="n", merge=T )

graphics.off()


#--------------------------------
#--------------------------------
#--------------------------------

climdat <- nc_open("/home/sscheiter/NcInputFiles/global_range/climvars_10min.nc")

xx <- ncvar_get( climdat, "lon" )
yy <- ncvar_get( climdat, "lat" )


prec <- NULL

for ( i in 1:length(x_coo) )
{
    for ( j in 1:length(y_coo) )
    {
        xmin <- which.min(abs(x_coo[i]-xx))
        ymin <- which.min(abs(y_coo[j]-yy))

        #print( c( x_coo[i], xx[xmin], x_coo[i]-xx[xmin] ) )
        #print( c( y_coo[j], yy[ymin], y_coo[j]-yy[ymin] ) )
	
        prec <- c( prec, sum( ncvar_get( climdat, "pre", c( xmin, ymin, 1 ), c( 1, 1, 12 ) ) ) )
    }
}



graphics.off()
pdf( "transect_map.pdf", height=8, width=14 )
par( mfcol=c(2,2), mar=c(4,4,1,1) )

plot( -100, -100, xlim=range(prec), ylim=c(0,10), xlab="MAP", ylab="Mean stem number" )
points( prec, c(t(sn0)), col="steelblue", pch=15 )
points( prec, c(t(sn1)), col="darkgoldenrod", pch=15 )

#--------------------------------

plot( -100, -100, xlim=range(prec), ylim=c(0,1), xlab="MAP", ylab="Tree:grass ratio" )
points( prec, c(t(vt0)), col="steelblue", pch=15 )
points( prec, c(t(vt1)), col="darkgoldenrod", pch=15 )

#--------------------------------

plot( -100, -100, xlim=range(prec), ylim=c(0,620), xlab="MAP", ylab="Tree woody biomass" )
points( prec, c(t(bw0)), col="steelblue", pch=15 )
points( prec, c(t(bw1)), col="darkgoldenrod", pch=15 )

#--------------------------------

plot( -100, -100, xlim=range(prec), ylim=c(0,2.5), xlab="MAP", ylab="Grass biomass" )
points( prec, c(t(bl0)), col="steelblue", pch=15 )
points( prec, c(t(bl1)), col="darkgoldenrod", pch=15 )

legend("topleft", c("No fire", "With fire"), col=c("steelblue", "darkgoldenrod"), lwd=2, bty="n", merge=T )

graphics.off()






