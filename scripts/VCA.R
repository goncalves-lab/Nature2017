###############################################################################
# Download the VCA assay data files (vca_*_dat.RData) from EGA.
# Source this file in R.
# See main().
###############################################################################

lmm_vars <- function( dat, facts, dir, s, e, update=F ) {
	library("lmerTest")
	
	dir.create( dir )
	nfile = paste( dir, "vars_", s, "-", e, ".RData", sep="")
	nfile2 = paste( dir, "fits_", s, "-", e, ".RData", sep="")
	
	if( !file.exists(nfile) & !update ) {
		vars = c()
		conv = c()
		fits = list()
		cg = NA
		for( i in s:e ) {
			cat( i, "\n" )
			ndat = cbind( y=dat[i,], facts )
			fits[[paste(i,"")]] = try( lmer( formula(paste( "y ~", paste( "( 1 |", colnames(facts), ")", collapse=" + " ) )), data=ndat ) )
			if( class(fits[[paste(i,"")]]) == "try-error" ) {
				vars = rbind(vars,rep(NA,ncol(facts)+1))
				tconv = NA
			} else {
				cg = i
				vars = rbind(vars,as.data.frame(VarCorr(fits[[paste(i,"")]]))$vcov)
				tconv = fits[[paste(i,"")]]@optinfo$conv$lme4$code
				if( is.null(tconv) ) {
					tconv = 0
				}
			}
			if( length(tconv) != 1 ) {
				cat( "\terror: tconv is", tconv, "\n" )
			}
			conv = c(conv,tconv[1])
		}
		colnames(vars) = as.data.frame(VarCorr(fits[[paste(cg,"")]]))$grp
		rownames(vars) = rownames(dat)[s:e]
		names(conv) = rownames(vars)
		cat( "saving variances and convergence to", nfile, "\n" )
		save( vars, conv, file=nfile )
		cat( "saving individual lmer fits to", nfile2, "\n" )
		save( fits, file=nfile2 )
	} else {
		cat( "file exists", nfile, "\n" )
	}
}

# join together the chunk results into a single file
lmm_vars_result <- function( dir ) {
	library( "svMisc" )
	library( "lmerTest" )
	files = dir( dir, paste("vars_[0-9]*-[0-9]*.RData", sep=""), full.names=T )
	
	m = max(as.numeric(sapply(strsplit(sapply(strsplit( sapply(strsplit(files,"vars_"),"[[",2), ".", fixed=T ), "[", 1), "-"), "[",2)))
	s = as.numeric(sapply(strsplit(sapply(strsplit(files[1],"vars_"),"[[",2),"-"),"[[",1))
	e = as.numeric(sapply(strsplit(sapply(strsplit(sapply(strsplit(files[1],"vars_"),"[[",2),"-"),"[[",2),".",fixed=T),"[[",1))
	
	avars = c()
	aconvs = c()
	cat( "reading files in", dir, "\n" )
	for( i in seq(1,m,(e-s+1)) ) {
		file = paste( dir, "vars_", i, "-", min(i+(e-s),m), ".RData", sep="" )
		if( file.exists(file) ) {
			cat( "\tloading", file, "\n" )
			load(file)
			avars = rbind(avars,vars)
			if( nrow(vars) != length(conv) ) {
				cat( "\t\terror: convergence vector has wrong length", length(conv), "should be", nrow(vars), "\n" )
			}
			aconvs = c( aconvs,conv ) 
		} else {
			cat( "\tfile does not exist", file, "\n" )
			empty = data.frame( row.names=seq(s,e) )
			for( j in 1:ncol(avars) ) {
				empty = cbind( empty, NA )
			}
			colnames(empty) = colnames(avars)
			avars = rbind(avars,empty)
			aconvs = c( aconvs, rep(NA,length(seq(s,e))) ) 
		}
	}	
	nfile = paste( dir, "vars.RData", sep="" )
	cfile = paste( dir, "conv.RData", sep="" )
	cat( "saving to", nfile, "\n" )
	save( avars, file=nfile )	
	save( aconvs, file=cfile )
}

run_vca <- function( assay ) {
	dir.create( assay )
	
	filename = paste("vca_", assay, "_dat.RData", sep="" )
	cat( "loading data from", filename, "\n" )
	load( filename )
	
	outdir = paste( assay, "/lmm_n", length(unique(facts$id_donor)), "/", sep="")
	
	# split into chunks that can be executed in parallel to speed up things
	ss = seq( 1, nrow(dat), e)
	for(s in ss) {
		lmm_vars( dat, facts, outdir, s, min(s+e-1,nrow(dat)) )
	}
	
	# join together the chunk results into a single file
	lmm_vars_result( outdir )
	
	# the results are in the variables 'avars' and 'aconvs'
	load(paste(outdir,"/vars.RData",sep=""))
	
	vpc = avars/rowSums(avars)
	return(vpc)
}

main <- function() {
	library("curl")
	
	setwd("your_directory_containing_the_downloaded_data_files/")
	
	url = "ftp://"
	qc1samples = read.table(curl(url,"r"),StringsAsFactors=F)[,1]
	url = "ftp://"
	qc2samples = read.table(curl(url,"r"),StringsAsFactors=F)[,1]
	
	vars = list()
	for( assay in c("gex","rnaseq","methyl","prot") ) {
		vars[[assay]] = run_vca("gex")
	}
}


