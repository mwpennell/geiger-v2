gbresolve=function(x, ...){
	UseMethod("gbresolve")
}

gbresolve.default=function(x, rank="phylum", within="", update=FALSE, split=FALSE){
	
	if(split) {
		tt=unique(tips<-sapply(x, function(y) unlist(strsplit(y, "_"))[1]))
	} else {
		tt=unique(tips<-x)
		names(tips)=tips
	}
	
	gb=.build.gbtaxdump(update)
	FUN=.fetch_gbhierarchy.above(gb, rank=rank, within=within)
	
	env=Sys.getenv()
	if("ignoreMULTICORE"%in%names(env)) {
		f=lapply
	} else {
		if(.check.multicore()) {
			f=function(X,FUN) mclapply(X,FUN,mc.silent=TRUE)
		} else {
			f=lapply
		}
	}
	
	tmp=f(tt, FUN)
	
	dd=sapply(tmp, function(y) all(is.na(y)))
	if(any(dd)) {
		warning(paste("The following taxa were not encountered in the NCBI taxonomy:\n\t", paste(tt[which(dd)], collapse="\n\t"), sep=""))
		tmp=tmp[-which(dd)]
		tt=tt[-which(dd)]
	}
	
	if(!length(tmp)) return(NULL)
	
	tmp=.compile_taxonomy(tmp)
	rownames(tmp)=tt
	
	res=matrix(tmp[match(tips, rownames(tmp)),], nrow=length(tips))
	rownames(res)=names(tips)
	colnames(res)=colnames(tmp)
	aa=apply(res, 1, function(x) all(is.na(x)))
	if(any(aa)) res=res[-which(aa),]
	res
}

## Assign internal node labels to phy based on genbank taxonomy
gbresolve.phylo=function(x, rank="phylum", within="", update=FALSE, split=TRUE){
	
	phy=x
	x=x$tip.label
	res=gbresolve(x, rank=rank, within=within, update=update, split=split)
	if(is.null(res)) return(res)
	ll=apply(res, 2, function(x) length(unique(x))==1 & !any(x==""))
	if(any(ll)){
		tmp=res[,1:min(which(ll))]
	} else {
		tmp=res
	}
	
	phy=nodelabel.phylo(phy, tmp)
	return(list(phy=phy, tax=res))
}

exemplar=function(x, ...) UseMethod("exemplar")

exemplar.phylo=function(x, strict.vals=NULL, taxonomy=NULL, ...){
	phy=x
# strict.vals: used if 
	if(is.null(taxonomy)) {
		tmp=gbresolve.phylo(phy, ...)
		taxonomy=tmp$tax
	}
	tax=taxonomy
	drp=which(!phy$tip.label%in%rownames(tax))
	z=exemplar(tax, strict.vals=strict.vals)
	ff=cbind(z, rownames(tax))
	ww=which(phy$tip.label%in%rownames(tax))
	phy$tip.label[ww]=z[match(phy$tip.label[ww], rownames(tax))]
	if(length(drp)) {
		phy=drop.tip(phy, phy$tip.label[drp])
	}
	phy
}

## FINDS MOST REPRESENTATIVE LINEAGE (i.e., exemplars) for EACH TAXON from a TAXONOMIC TABLE
#Chioglossa_lusitanica           "Chioglossa"       "Salamandridae"   -->	Chioglossa
#Lyciasalamandra_antalyana       "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_antalyana   
#Lyciasalamandra_atifi           "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_atifi       
#Lyciasalamandra_billae          "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_billae      
#Lyciasalamandra_fazilae         "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_fazilae     
#Lyciasalamandra_flavimembris    "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_flavimembris
#Lyciasalamandra_helverseni      "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_helverseni  
#Lyciasalamandra_luschani        "Lyciasalamandra"  "Salamandridae"   -->	Lyciasalamandra_luschani    
#Mertensiella_caucasica          "Mertensiella"     "Salamandridae"   -->	Mertensiella      
#Salamandra_algira               "Salamandra"       "Salamandridae"   -->	Salamandra_algira           
#Salamandra_atra                 "Salamandra"       "Salamandridae"   -->	Salamandra_atra             
exemplar.default=function(x, strict.vals=NULL){
	
	tax=x
	if(!is.matrix(tax) | is.null(rownames(tax)) | is.null(colnames(tax))) stop("supply 'tax' as a matrix with both row and column names")
	if(!is.null(strict.vals)){
		if(!all(is.character(strict.vals))) stop("supply 'strict.vals' as a character vector specifying which columns are to be used in 'tax'")
		nn=sort(match(strict.vals, colnames(tax)))
		nn=nn[!is.na(nn)]
		if(length(nn)){
			tax=as.matrix(tax[,nn])
		}
	}
	incomparables=c("", NA)
	z=rownames(tax)
	for(j in 1:ncol(tax)){
		tt=table(tax[,j])	
		tt=tt[!names(tt)%in%incomparables]
		if(any(jj<-(tt==1))){
			nn=names(tt[jj])
			mm=match(nn, tax[,j])
			z[mm]=nn
		}	
	}
	z
}

.fetch_gbhierarchy.above=function(gb, rank="root", within=""){
	
## returns taxonomic information for a 'taxon' up to the 'rank' given
## requires fetch_genbank.pl and potentially nodes.dmp and names.dmp (if /tmp/idx is not available)
	
	Linnaean=c(
		       "root",
			   "superkingdom",
			   "kingdom",
			   "superphylum",
			   "phylum",
			   "subphylum",
			   "superclass",
			   "class",
			   "subclass",
			   "infraclass",
			   "superorder",
			   "order",
			   "suborder",
			   "parvorder",
			   "infraorder",
			   "superfamily",
			   "family",
			   "subfamily",
			   "tribe",
			   "subtribe",
			   "genus",
			   "subgenus",
#			   "species group",
#		       "species subgroup",
			   "species",
			   "subspecies",
			   "varietas",
			   "forma"
	)
	
	dat=gb	
	
	get_tax=function(name){
		
		# resolve highest rank requested
		if(!is.null(rank)) if(!rank%in%Linnaean) {
			cat(paste("Ensure that supplied 'rank' is in: \n\t", paste(Linnaean, collapse="\n\t"), "\n", sep=""))
			stop("Supplied 'rank' is unrecognized.")
		}
		
		rank=match.arg(rank, Linnaean)
		
		name=gsub("_", " ", name) 
		sci=dat$type=="scientific name"

		fetch_anc=function(id){
			ww=which(dat$id==id & sci==TRUE)
			if(length(ww)!=1) return(NA)
			ww
		}
		

		# resolve 'name' to 'id'
		idx=which(dat$node==name)
		
		## NO MATCH -- use grep
		if(!length(idx)){
			idx=grep(name, dat$node, ignore.case=TRUE)
			if(!length(idx)==1){
				if(length(unique(dat$parent_id))==1){
					idx=min(idx)
				} else {
					if(length(idx)>1) warning(paste("Attempt one of the following:\n\t", paste(dat$node[idx], collapse="\n\t"), sep=""))
					return(NA)
				}
			}
		}
		
		## MORE THAN ONE MATCH -- see if all point to same 'anc'
		if(length(idx)>1){
			if(length(unique(dat$parent_id))==1){
				idx=min(idx)
			} else {
				if(within=="") {
					warning(paste(sQuote(name), "does not appear to be a unique 'name'"))
					return(NA)
				}
			}
		}
		
		alltax=function(idx){
			## COMMON NAME --- resolve to scientific name
			if(!all(dat$type[idx]=="scientific name")){
				idx=sapply(idx, function(x){
					while(1){
						   xx=fetch_anc(dat$id[x])
						   if(dat$type[xx]=="scientific name") break()
					}
					return(xx)
				})
				
				
			}
			
			id=dat$id[idx]
			
			
			## COLLECT TAXONOMY
			pnms=character()
			rnks=character()
			
			cur=id
			rr=""
			while(rr!=rank){
				orig=cur
				idx=fetch_anc(cur)
				pnm=dat$node[idx]
				pid=dat$parent_id[idx]
				pnms=c(pnms, pnm)
				cur=pid
				rr=dat$rank[idx]
				rnks=c(rnks, rr)
				if(rr==rank || orig==pid || pnm=="root") break()
			}
			
			res=pnms
			names(res)=rnks
			res=res[names(res)%in%Linnaean]
			return(res)			
		}
		
		if(length(idx)>1){
			tmp=lapply(idx, alltax)
			ww=sapply(tmp, function(x) tolower(within)%in%tolower(x))
			if(any(ww)){
				if(sum(ww)>1){
					uu=table(names(unlist(tmp[which(ww)])))
					uu=names(uu[uu==sum(ww)])
					x=sapply(tmp[which(ww)], function(y) digest(y[names(y)%in%uu]))
					if(length(unique(x))==1){
						ll=sapply(tmp[which(ww)], length)
						return(tmp[[which(ww)[min(which(ll==max(ll)))]]])
					} else {
						warning(paste("Attempt one of the following:\n\t", paste(sapply(tmp[which(ww)], function(x) x[[1]]), collapse="\n\t"), sep=""))
					}
				} else {
					return(tmp[[which(ww)]])
				}
			} else {
				return(NA)
			}
		} else {
			if(within!=""){
				tmp=alltax(idx)
				if(tolower(within)%in%tolower(tmp)){
					return(tmp)
				} else {
					warning(paste(sQuote(name), "was encountered in NCBI taxonomy but not found within", sQuote(within)))
					return(NA)

				}
			} else {
				tmp=alltax(idx)
				if(!rank%in%names(tmp) & rank!="root"){
					warning(paste(sQuote(rank), "appears to be a rank inconsistent with NCBI taxonomy for", sQuote(name)))
				}
				return(tmp)
			}
		}
	}
	
	return(get_tax)
}

.path_to_gb.dmp=function(package="geiger"){
	path=paste(system.file(package=package), "data", sep="/")
#	pl=paths[[which(names(paths)=="fetch_genbank")]]
#	base_path=gsub("fetch_genbank.pl", "", pl)
	nd=paste(path, "nodes.dmp", sep="/")
	nm=paste(path, "names.dmp", sep="/")
	ncbi=paste(path, "_ncbi.rda", sep="/")
	return(c(ncbi=ncbi, names=nm, nodes=nd, base=path))
}


## EXPORTING GB taxonomy to R tables
# partly from OMeara phyloorchard code ncbiTaxonomy.R::ncbiTaxonomy()
.build.gbtaxdump=function(update=FALSE){
	
	gb_path=.path_to_gb.dmp()
	rda=as.list(gb_path)$ncbi
	
	build=update
	if(!file.exists(rda) | build){
		build=TRUE
	} else {
		build=FALSE
	}

	if(build){
		if(file.exists("taxdump.tar.gz")) unlink("taxdump.tar.gz")
		cat("Please be patient as 'taxdump' is built from NCBI; download may take several minutes...\n")
		if(!system("which wget", ignore.stdout=TRUE)==0) stop("Install 'wget' before proceeding.")
		if(!system("which gunzip", ignore.stdout=TRUE)==0) stop("Install 'gunzip' before proceeding.")
		if(!system("which tar", ignore.stdout=TRUE)==0) stop("Install 'tar' before proceeding.")
		if(!system("which perl", ignore.stdout=TRUE)==0) stop("Install 'perl' before proceeding.")
		
		
		system("wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
		system("gunzip taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
		system("tar -xvf taxdump.tar", ignore.stderr=TRUE, wait=TRUE)
		system("perl -i -p -e's/[^\\w+^\\s+^\\d+^\\|+]//g' names.dmp")
		system("perl -i -p -e's/[^\\w+^\\s+^\\d+^\\|+]//g' nodes.dmp")
		
		mv=TRUE
		rm=TRUE
		
		if(mv){
			system("rm -rf /tmp/idx")
			if(all(sapply(ff<-c("nodes.dmp", "names.dmp"), file.exists))){
				system(paste("mv nodes.dmp", gb_path[["nodes"]], sep=" "))
				system(paste("mv names.dmp", gb_path[["names"]], sep=" "))
				cleanup=c("taxdump.tar", "readme.txt", "gc.prt", "merged.dmp", "division.dmp", "delnodes.dmp", "citations.dmp", "gencode.dmp")
				cleanup=cleanup[cleanup%in%dir()]
				if(rm) system(paste("rm -f", paste(cleanup, collapse=" "), sep=" "))
			} else {
				stop("Error encountered from NCBI: 'nodes.dmp' and (or) 'names.dmp' cannot be located.")
			}
			
			## CALL UP gb data
			names.dmp<-read.table(gb_path[["names"]],header=FALSE, sep="|",strip.white=TRUE,fill=TRUE,stringsAsFactors=FALSE) 
			names.dmp<-names.dmp[,1:4]
			names(names.dmp)<-c("id", "node", "unique", "type")
			nodes.dmp<-read.table(gb_path[["nodes"]],header=FALSE, sep="|",strip.white=TRUE,fill=TRUE,stringsAsFactors=FALSE)
			nodes.dmp<-nodes.dmp[,c(1:5,11:13)]
			names(nodes.dmp)<-c("id","parent_id","rank","embl_code","division_id","GenBank_hidden_flag","hidden_subtree_root_flag","comments")

#		gb=list(nodes=nodes.dmp, names=names.dmp)
#			class(gb)=c("taxdump", class(gb))
#			save(gb, file=rda)
			ncbi=cbind(names.dmp, nodes.dmp[match(names.dmp$id, nodes.dmp$id), c("parent_id", "rank")])
			rownames(ncbi)=NULL
			attr(ncbi, "date")=Sys.Date()
			class(ncbi)=c("taxdump", class(ncbi))
			save(ncbi, file=rda)
			
		}
	}
	

	ncbi=get(load(rda))
	return(ncbi)

}

#require(auteur)
#data(urodela)
#phy=urodela$phy


## Assign internal node labels to phy based on genbank taxonomy
gbresolve.phylo=function(x, rank="phylum", within="", update=FALSE, split=TRUE){
	phy=x
#	require(auteur)
	if(split) {
		tt=unique(tips<-sapply(phy$tip.label, function(x) unlist(strsplit(x, "_"))[1]))
	} else {
		tt=unique(tips<-phy$tip.label)
		names(tips)=tips
	}
	
	gb=.build.gbtaxdump(update)
	FUN=.fetch_gbhierarchy.above(gb, rank=rank, within=within)
	
	env=Sys.getenv()
	if("ignoreMULTICORE"%in%names(env)) {
		f=lapply
	} else {
		if(.check.multicore()) {
			f=function(X,FUN) mclapply(X,FUN,mc.silent=TRUE)
		} else {
			f=lapply
		}
	}
	tax=f(tt, function(x) try(FUN(x), silent=TRUE))
	dd=sapply(tax, function(x) inherits(x, "try-error") | all(is.na(x)))
	
	if(any(dd)){
		tt=tt[!dd]
		tax=tax[!dd]
	}
	
	names(tax)=tt
	
	if(!length(tt)) return(NA)
	
	tmp=.compile_taxonomy(tax)
	res=tmp[match(tips, rownames(tmp)),]
	rownames(res)=names(tips)
	aa=apply(res, 1, function(x) all(is.na(x)))
	if(any(aa)) res=res[-which(aa),]
	ll=apply(res, 2, function(x) length(unique(x))==1 & !any(x==""))
	if(any(ll)){
		tmp=res[,1:min(which(ll))]
	} else {
		tmp=res
	}
	
	phy=nodelabel.phylo(phy, tmp)
	return(list(phy=phy, tax=res))
}

print.taxdump=function(x, ...){
	cat(paste("\nNCBI GenBank taxonomy assembled ", attributes(x)$date, "\n ...showing the first several entries...\n\n", sep=""))	
	print(head(as.data.frame(x)))
	
}

.compile_taxonomy=function(tax){
	all=c("forma", "superkingdom", "kingdom", 
		  "superphylum", "phylum", "subphylum", "superclass", "class", 
		  "subclass", "infraclass", "superorder", "order", "suborder", 
		  "parvorder", "infraorder", "superfamily", "family", "subfamily", 
		  "tribe", "subtribe", "genus", "subgenus", "species", 
		  "subspecies", "varietas")
	tmp=sapply(tax, names)
	drop=sapply(tmp, function(x) all(is.null(x)))
	dat=tax
	if(any(drop)) {
		dat=dat[-which(drop)]
		tmp=tmp[-which(drop)]
	}
	taxa=tmp
	levels=unique(unlist(taxa))
	mm=match(all, levels)
	hier=rev(all[!is.na(mm)])
	mm=matrix("", nrow=length(dat), ncol=length(hier))
	for(i in 1:length(dat)){
		cur=dat[[i]]
		mm[i,match(names(cur), hier)]=cur
	}
	
	## exclude peculiar assignments (e.g., Proteus [salamander genus] returns prokaryotic labels)
	## FIX ME: more problem anticipation
#	primary=table(mm[,ncol(mm)])
#	zz=mm[,ncol(mm)]!=names(primary[primary==max(primary)])
#	if(any(zz)) {
#		mm=mm[-which(zz),]
#		dat=dat[-which(zz)]
#	}
	mm[is.na(mm)]=""
	rownames(mm)=names(dat)
	colnames(mm)=hier
	mm
}


