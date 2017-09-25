#############################
## Load library / function ##
#############################

computeGenotypes <- Vectorize( function( ind, snp, Prob ) {
	prob_AA <- Prob[ ind, snp ]
	prob_AB <- Prob[ ind + 1, snp ]
	prob_BB <- Prob[ ind + 2, snp ]
	prob_sum <- prob_AA + prob_AB + prob_BB
	( prob_AB + 2 * prob_BB ) / prob_sum
}, 'ind' )

####################
## Read GEN file  ##
## Assuming there ##
## is chromo col  ##
####################

file <- 'genfile.gen'

lines <- strsplit( readLines( file, n = 1 ), split = ' ' )[[ 1 ]]
what <- vector( 'list', length( lines ) )
what[[ 1 ]] <- integer( 0 )
what[[ 2 ]] <- character( 0 )
what[[ 3 ]] <- character( 0 )
what[[ 4 ]] <- integer( 0 )
what[[ 5 ]] <- character( 0 )
what[[ 6 ]] <- character( 0 )
for( ii in 7 : length( lines ) ) {
	what[[ ii ]] <- numeric( 0 )
}

res <- scan( file, what = what )
Meta <- cbind( as.data.frame( res[ seq_len( 6 ) ], stringsAsFactors = FALSE ) )
colnames( Meta ) <- c( 'chromosome', 'SNPID', 'RSID', 'position', 'A_allele', 'B_allele' )
Prob <- do.call( 'rbind', res[ -seq_len( 6 ) ] )

################################
## Compute correlation matrix ##
################################

N <- nrow( Prob ) / 3
m <- ncol( Prob )

G <- t( computeGenotypes( seq( 1, 3 * N, 3 ), seq_len( m ), Prob ) )
R <- cor( G )

##################
## Write output ##
##################

write.table( R, file = sub( '.gen', '.ld', file ), row.names = F, col.names = F, quote = F )