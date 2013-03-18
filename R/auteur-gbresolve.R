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
#   "species group",
#   "species subgroup",
	"species",
	"subspecies",
	"varietas",
	"forma"
)



.path_to_gb.dmp=function(package="geiger"){
	path=paste(system.file(package=package), "data", sep="/")
#	pl=paths[[which(names(paths)=="fetch_genbank")]]
#	base_path=gsub("fetch_genbank.pl", "", pl)
	nd=paste(path, "nodes.dmp", sep="/")
	nm=paste(path, "names.dmp", sep="/")
	ncbi=paste(path, "ncbi.rda", sep="/")
	return(c(ncbi=ncbi, names=nm, nodes=nd, base=path))
}


## EXPORTING GB taxonomy to R tables
# partly from OMeara phyloorchard code ncbiTaxonomy.R::ncbiTaxonomy()
.ncbi=function(update=FALSE){
	
	gb_path=.path_to_gb.dmp()
	rda=as.list(gb_path)$ncbi
	
	build=update
	if(!file.exists(rda) | build){
		build=TRUE
	} else {
		build=FALSE
		if(exists("gbtaxdump")) return(gbtaxdump)
	}
	
	if(build){
		if(file.exists("taxdump.tar.gz")) unlink("taxdump.tar.gz")
		cat("Please be patient as 'taxdump' is built from NCBI; download may take several minutes...\n")
		if(!system("which curl", ignore.stdout=TRUE)==0) stop("Install 'curl' before proceeding.")
		if(!system("which gunzip", ignore.stdout=TRUE)==0) stop("Install 'gunzip' before proceeding.")
		if(!system("which tar", ignore.stdout=TRUE)==0) stop("Install 'tar' before proceeding.")
		if(!system("which perl", ignore.stdout=TRUE)==0) stop("Install 'perl' before proceeding.")
		
		
		system("curl -OL ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
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
            ncbi[which(ncbi[,"node"]=="root"),"rank"]="root"
			attr(ncbi, "date")=Sys.Date()
			class(ncbi)=c("taxdump", class(ncbi))
			save(ncbi, file=rda)
			
		}
	}
	
	
	ncbi=get(load(rda))
	return(ncbi)
	
}

gbtaxdump=.ncbi(update=FALSE)

gbcontain=function(x, rank="species", within="", ...){
	type="scientific name"
	rank=match.arg(rank, Linnaean)
	lidx=which(Linnaean==rank)
    if(is.na(lidx)) stop("uninterpretable 'rank'")
#	if(length(x)>1) stop("Supply 'x' as a single taxon")
	
	gb=.ncbi(...)
	ss=gb[,"type"]%in%type
#    gb=cbind(gb, ridx=rr)
    gb=gb[ss,]
#    gb=gb[gb$ridx<=lidx,]
    

	#FUN=function(id){
	#	keep=c()
    #		wx=which(gb[,"parent_id"]==id & ss)
	#	if(length(wx)){
	#		keep=c(keep, gb[wx[drop<-which(gb[wx,"rank"]==rank)],"node"])
	#		if(length(drop)) wx=wx[-drop]
	#		if(length(wx)) keep=c(keep, unlist(lapply(gb[wx,"id"], FUN)))
	#	}
	#	return(keep)
	#}
    
    FUN=function(id){
        rr=match(gb[,"rank"], Linnaean)

        k=rep(1,nrow(gb))
        k[which(rr==lidx)]=0
        dfind=.gb_des_worker(gb$parent_id, gb$id, k)
        ids=dfind(id)
        ww=match(ids, gb$id)
        nm=gb[ww[gb$rank[ww]%in%rank],"node"]
        nm
    }
	
	gbc=function(x){
		if(is.character(x)) {
			ww=which(tolower(gb[,"node"])==tolower(x))
		} else if(is.numeric(x)){
			ww=which(gb[,"id"]==x)
		} else {
			return(NULL)
		}
		
		if(length(ww)==1){
			id=gb[ww,"id"]
			return(FUN(id))
		} else if(!length(ww)){
			return(NULL)
		} else {
			if(within!=""){
				par=gb[ww, "parent_id"]
				tmp=gbresolve(par)
				if(!"matrix"%in%class(tmp)) return(NULL)
				gg=apply(tmp, 1, function(x) any(x%in%within))
				if(sum(gg)==1){
					id=gb[ww[which(gg)],"id"]
					return(gbc(id))
				} else {
					return(NULL)
				}
			}
			return(x)
		}
	}
	
	
	env=Sys.getenv()
	if("ignoreMULTICORE"%in%names(env)) {
		f=lapply
	} else {
		f=.get.multicore()
	}
	
	res=f(x, function(g) {gbc(g)})

	probs=sapply(1:length(res), function(idx) all(x[idx]==res[[idx]]))
	
	if(any(probs)){
		warning(paste("Try using the 'within' argument as the following taxa are not unique:\n\t", paste(x[probs], collapse="\n\t"), sep=""))
	}
	
	res=res[!probs]
	if(length(res)){
		names(res)=x[!probs]
		return(res)
	} else {
		return(NULL)
	}
}


.gb_des_worker=function(anc, des, keep){
    fetcher=function(node){
        if(length(unique(c(length(anc), length(des)))->sz)!=1) stop("'anc' and 'des' appear mismatched")
		mat=list(
            node=as.integer(node),
            nrow=as.integer(sz),
            ANC=as.integer(anc),
            DES=as.integer(des),
            keep=as.integer(keep)
        )
		.Call("compile_descendants", mat=mat, package="geiger")[[1]]
    }
    fetcher
}


gbresolve=function(x, rank="phylum", within="", ...){
	UseMethod("gbresolve")
}

gbresolve.default=function(x, rank="phylum", within="", ...){
		
    ridx=match(rank, Linnaean)
    if(any(is.na(ridx)) | length(ridx)>2) stop("'rank' should be a vector of one or two elements occuring in Linnaean") 
    rank=rank[oridx<-order(ridx)]
    ridx=ridx[oridx]

	gb=.ncbi(...)

	FUN=.fetch_gbhierarchy.above(gb, rank=rank[1], within=within)
				
	if(all(is.numeric(x) | is.integer(x))){
		ss=gb[,"type"]=="scientific name"
		x=sapply(x, function(y) gb[which(gb[,"id"]==y & ss),"node"])				  
	}
		

    tt=unique(tips<-x)
    names(tips)=tips
						  
	
	f=.get.multicore()
	
	tmp=f(tt, FUN)
	
	dd=sapply(tmp, function(y) all(is.na(y)))
	if(any(dd)) {
		warning(paste("The following taxa were not encountered in the NCBI taxonomy:\n\t", paste(tt[which(dd)], collapse="\n\t"), sep=""))
		tmp=tmp[-which(dd)]
		tt=tt[-which(dd)]
	}
    
    tmp=lapply(1:length(tmp), function(idx) {
        x=tmp[[idx]]
        mm=match(names(x), Linnaean)
        if(length(ridx)==2) ss = mm>=ridx[1] & mm<=ridx[2] else ss = mm>=ridx[1]
        if(!any(ss)) return(NA)
        y=x[ss]
        y
    })
    
    if(!length(tmp)) return(NULL)

    dd=sapply(tmp, function(y) all(is.na(y)))
    if(any(dd)) {
		warning(paste("The following taxa were not encountered in the NCBI taxonomy:\n\t", paste(tt[which(dd)], collapse="\n\t"), sep=""))
		tmp=tmp[-which(dd)]
		tt=tt[-which(dd)]
	}
	
	if(!length(tmp)) return(NULL)
	
	names(tmp)=tt
	tmp=.compile_taxonomy(tmp)
	
	res=matrix(tmp[match(tips, rownames(tmp)),], nrow=length(tips))
	rownames(res)=names(tips)
	colnames(res)=colnames(tmp)
	aa=apply(res, 1, function(x) all(is.na(x)))
	if(any(aa)) res=res[-which(aa),]
	res
}

## Assign internal node labels to phy based on genbank taxonomy
gbresolve.phylo=function(x, rank="phylum", within="", ...){
	
	phy=x
	x=x$tip.label
	res=gbresolve(x, rank=rank, within=within, ...)
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

subset.phylo=function(x, taxonomy, rank="", ...){
## rank (e.g., 'family') and 'family' must be in columns of 'taxonomy'
	
	phy=x
	if(!rank%in%colnames(taxonomy)){
		stop(paste(sQuote(rank), " does not appear as a column name in 'taxonomy'", sep=""))
	}
	
	xx=match(phy$tip.label, rownames(taxonomy))
	
	new=as.matrix(cbind(tip=phy$tip.label, rank=taxonomy[xx,rank]))
	drop=apply(new, 1, function(x) if( any(is.na(x)) | any(x=="")) return(TRUE) else return(FALSE))
	if(any(drop)){
		warning(paste("Information for some tips is missing from 'taxonomy'; offending tips will be left unpruned:\n\t", paste(phy$tip.label[drop], collapse="\n\t"), sep=""))
#		phy=drop.tip(phy, phy$tip.label[drop])
#		new=new[!drop,]
	}
	
	tips=phy$tip.label
	hphy=hashes.phylo(phy, tips=tips)
	tax=as.data.frame(new, stringsAsFactors=FALSE)
	stax=split(tax$tip,tax$rank)
	rank_hashes=sapply(stax, function(ss) .hash.tip(ss, tips=tips))
	
	pruned=hphy
	pruned$tip.label=ifelse(drop==TRUE, tax$tip, tax$rank)

	if(!all(zz<-rank_hashes%in%hphy$hash)){
		warning(paste(paste("non-monophyletic at level of ",rank,sep=""),":\n\t", paste(sort(nonmon<-names(rank_hashes)[!zz]), collapse="\n\t"), sep=""))
#       for(j in 1:length(nonmon)){
            vv=which(pruned$tip.label%in%nonmon)
            pruned$tip.label[vv]=tax$tip[vv]
#       }
#		pruned=drop.tip(pruned, nonmon)
        
	}
		
	rank_phy=unique.phylo(pruned)
	rank_phy$tip.label=as.character(rank_phy$tip.label)
	return(rank_phy)	
}


exemplar.phylo=function(phy, taxonomy=NULL, ...){
	if(is.null(taxonomy)) {
		tmp=gbresolve.phylo(phy, ...)
		taxonomy=tmp$tax
	}
	tax=taxonomy
	drp=which(!phy$tip.label%in%rownames(tax))
	z=.exemplar.default(tax)
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
.exemplar.default=function(x){
	
	tax=x
	if(!is.matrix(tax) | is.null(rownames(tax)) | is.null(colnames(tax))) stop("supply 'tax' as a matrix with both row and column names")

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

.gb_anc_worker=function(anc, des, root){
	fetcher=function(node){
		if(length(unique(c(length(anc), length(des)))->sz)!=1) stop("'anc' and 'des' appear mismatched")
		mat=list(
            node=as.integer(node),
            root=as.integer(root),
            nrow=as.integer(sz),
            ANC=as.integer(anc),
            DES=as.integer(des)
        )
		.Call("compile_ancestors", mat=mat, package="geiger")[[1]]
	}
	fetcher
}


### NEEDS TO BE WRITTEN TO TAKE ADVANTAGE OF .GB_ANC_WORKER
.fetch_gbhierarchy.above=function(gb, rank="root", within=""){
	
## returns taxonomic information for a 'taxon' up to the 'rank' given
	dat=gb	
    sci=dat$type=="scientific name"
    gb=gb[sci,]

    ancFUN=.gb_anc_worker(gb$parent_id, gb$id, 1)

	get_tax=function(name){
		
		# resolve highest rank requested
		if(!is.null(rank)) if(!rank%in%Linnaean) {
			cat(paste("Ensure that supplied 'rank' is in: \n\t", paste(Linnaean, collapse="\n\t"), "\n", sep=""))
			stop("Supplied 'rank' is unrecognized.")
		}
		
		rank=match.arg(rank, Linnaean)
		
		name=gsub("_", " ", name) 

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
#				if(within=="") {
#					warning(paste(sQuote(name), "does not appear to be a unique 'name'"))
#					return(NA)
#				}
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
			
            tmp=ancFUN(id)
            mm=match(tmp, gb$id)
            pnms=gb$node[mm]
            rnks=gb$rank[mm]
            res=pnms
            names(res)=rnks
    
			return(res[rnks%in%Linnaean])			
		}
		
		if(length(idx)>1){
			tmp=lapply(idx, alltax)
			ww=sapply(tmp, function(x) tolower(within)%in%tolower(x))
			if(any(ww) | within==""){
                if(within=="") ww=sapply(tmp, length)>0
				if(sum(ww)>1 | within==""){
					uu=table(names(unlist(tmp[which(ww)])))
					uu=names(uu[uu==sum(ww)])
					x=sapply(tmp[which(ww)], function(y) digest(y[names(y)%in%uu]))
					if(length(unique(x))==1){
						ll=sapply(tmp[which(ww)], length)
						return(tmp[[which(ww)[min(which(ll==max(ll)))]]])
					} else {
						warning(paste("Attempt one of the following:\n\t", paste(sapply(tmp[which(ww)], function(x) x[[1]]), collapse="\n\t"), sep=""))
                        return(NA)
					}
				} else {
					return(tmp[[which(ww)]])
				}
			} else {
                if(within!="") {
                    if(length(tmp)) warning(paste(sQuote(name), "was encountered in NCBI taxonomy but not found within", sQuote(within)))
                }
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


print.taxdump=function(x, ...){
	if(nrow(x)>6){
		cat(paste("\nNCBI GenBank taxonomy assembled ", attributes(x)$date, "\n ...showing the first several entries...\n\n", sep=""))	
	} else {
		cat(paste("\nNCBI GenBank taxonomy assembled ", attributes(x)$date, "\n\n", sep=""))	
	}
	print(head(as.data.frame(x)))
	
}

.compile_taxonomy=function(tax){
	all=Linnaean	
	
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
	
	mm[is.na(mm)]=""
	rownames(mm)=names(dat)
	colnames(mm)=hier
	mm
}


