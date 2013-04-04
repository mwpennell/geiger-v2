
## GENERIC FUNCTION
heights <- function(x)
UseMethod("heights")


heights.phylo=function(x){
    phy=x
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root = ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs = c(phy$tip.label, phy$node.label)
	depth = max(xx)
	tt = depth - xx
	idx = 1:length(tt)
	dd = phy$edge.length[idx]
	mm = match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd = c(phy$edge.length, root)[mm]
	ss = tt + dd
	res = cbind(ss, tt)
	rownames(res) = idx
	colnames(res) = c("start", "end")
	res = data.frame(res)
	res	
}


heights.multiPhylo=function(x){
    phy=x
	phy=hashes.phylo(phy)
	hh=unique(unlist(lapply(phy, function(x) x$hash)))
	times=lapply(phy, function(x) {tmp=heights.phylo(x); tmp$hash=x$hash; tmp})
	if(is.null(names(times))) names(times)=1:length(times)
	out=lapply(hh, function(hash){
			   tmp=sapply(times, function(dat){
						  if(hash%in%dat$hash) dat[which(dat$hash==hash),c("start","end")]  else return(c(0,0))
						  })
			   res=data.frame(t(tmp))
			   rownames(res)=names(times)
			   names(res)=c("start","end")
			   gg=apply(res, 1, function(x) all(x==0))
			   if(any(gg)) res=res[-which(gg),]
			   return(res)
			   })
	attr(out,"hash")=hh
	out
}

.unique.phylo=function(phy){
# phy is assumed to have tips that are redundant (and whose exemplars are monophyletic)
# prunes tree to leave one member of each unique label
	
	mon=table(phy$tip.label)
	if(all(mon==1)) return(phy)
	mon=mon[mon>1]
	for(i in 1:length(mon)){
		m=names(mon[i])
		drop=which(phy$tip.label==m)
		drop=drop[2:length(drop)]
		phy$tip.label[drop]="null"
	}
	phy=drop.tip(phy, phy$tip.label[phy$tip.label=="null"])
	return(phy)
}


unique.phylo=function(x, incomparables=FALSE, ...){
# phy: has tip.labels that are not unique
# returns phylogeny with unique tip labels
# if a taxon (indicated by multiple tip labels of same string) is non-monophyletic, all tips of that taxon are removed with warning
# likely used where an exemplar phylogeny is needed and phy$tip.label is 'tricked' into a non-unique vector of names
	
	phy=x
	if(incomparables) warning("'incomparables' exerts no effect in this setting")
	.exemplar.monophyly=function(tip, phy) {
# tip: occurs multiply in phy$tip.label
		nn=which(phy$tip.label==tip)
		if(length(nn)==1) return(TRUE) 
		anc=getMRCA(phy, nn)
		dd=unique(phy$tip.label[.get.descendants.of.node(anc, phy, tips=TRUE)])
		if(length(dd)==1) return(TRUE) else return(FALSE)
	}

	
	tt=table(phy$tip.label)
	if(!any(tt>1)) {
		return(phy)
	} else {
		todo=names(tt[tt>1])
#		cat("checking monophyly of groups...\n")
		mon=sapply(todo, function(t) {
#				   cat(paste("\n\t",t, sep=""))
				   .exemplar.monophyly(t, phy)
				   })
#		cat("\n")
		
		if(any(!mon)) {
			warning(paste("non-monophyletic lineages encountered:\n\t", paste(todo[!mon], collapse="\n\t"), "\n", sep=""))
			phy=drop.tip(phy, names(!mon))
		}
		return(.unique.phylo(phy))
	}
}

unique.multiPhylo=function(x, incomparables=FALSE, ...){
	phy=x
	if(incomparables) warning("'incomparables' exerts no effect in this setting")
	ss=sapply(phy, digest)
	if(any(dd<-duplicated(ss))){
		sub=phy[-which(dd)]
		class(sub)="multiPhylo"
		return(sub)
	} else {
		return(phy)
	}
}



#general phylogenetic utility for determining whether a node is the root of the phylogeny
#author: JM EASTMAN 2011

is.root <-
function(node,phy) {
	if(node==Ntip(phy)+1) return(TRUE) else return(FALSE)
}


.nodefind.phylo=function(phy, label){
	N=Ntip(phy)
	ww=match(label, phy$node.label)
	if(!all(is.na(ww))) {
		return(N+ww)
	} else {
		warning(paste("Encountered no node.label for ",label,sep=""))
		return(NULL)
		
	}
}

.cache.descendants=function(phy){
# fetches all tips subtended by each internal node
	
	N=as.integer(Ntip(phy))
	n=as.integer(Nnode(phy))
	
	phy=reorder(phy, "postorder")
	
	zz=list( N=N,
			MAXNODE=N+n,
			ANC=as.integer(phy$edge[,1]), 
			DES=as.integer(phy$edge[,2])
			)
	
	res=.Call("cache_descendants", phy=zz, package="geiger")
	return(res)
}


.vmat <- function(phy){
	n=Ntip(phy)
	out <- .Call("vmat", tree=list(
								   ROOT = as.integer(n+1),
								   MAXNODE = as.integer(max(phy$edge[,1])),
								   ENDOFCLADE = as.integer(dim(phy$edge)[1]),
								   ANC = as.integer(phy$edge[,1]),
								   DES = as.integer(phy$edge[,2]),
								   EDGES = as.double(c(phy$edge.length,0)),
								   VCV = as.double(array(matrix(0, n, n)))),
				 PACKAGE = "geiger")
	v=matrix(out$VCV,nrow=n,byrow=FALSE)
	rownames(v)<-colnames(v)<-phy$tip.label
	return(v)
} 


.polytomy.phylo=function(tips, age=1){
	N=length(tips)+1
	edge=cbind(N,1:(N-1))
	length=rep(age, N-1)
	phy=list(edge=edge, edge.length=length, tip.label=tips, Nnode=1)
	class(phy)="phylo"
	return(phy)
}


nodelabel.phylo=function(phy, taxonomy, strict=TRUE){
# all phy$tip.label must be in taxonomy
# taxonomy: exclusivity highest on left, lowest on right (species, genus, family, etc., as columns)
# columns in 'taxonomy' should ONLY be taxonomic ranks
	

#	tt=as.matrix(taxonomy)
#	tt[!is.na(tt)]=""
#	drp=apply(tt, 1, function(x) all(x==""))
#	if(any(drp)) taxonomy=taxonomy[-which(drp),]
	
	
	taxonomy=cbind(rownames=rownames(taxonomy),taxonomy)
	rank="rownames"
	
	taxonomy=as.data.frame(as.matrix(taxonomy),stringsAsFactors=FALSE)
			
	if(!all(xx<-phy$tip.label%in%taxonomy[,rank])) {
		warning(paste("taxa not found in 'taxonomy':\n\t", paste(phy$tip.label[!xx], collapse="\n\t"), sep=""))
	}
	taxonomy[taxonomy==""]=NA

    op<-orig<-options()
    op$expressions=max(op$expressions, 500000)
	options(op)
	
	unmatched=phy$tip.label[!xx]
	idx=match(phy$tip.label, taxonomy[,rank])
	tax=taxonomy[idx,]
	
	labels=unique(unlist(tax[,-which(names(tax)==rank)]))
	labels=labels[!is.na(labels)]
	
	dat=tax[,-which(names(tax)==rank)]
	hashes_labels=character(length(labels))
	zz=tax[,rank]
	tips=phy$tip.label[xx]
	for(i in 1:ncol(dat)){
		uu=unique(dat[,i])
		uu=uu[!is.na(uu)]
		for(j in uu){
			cur=zz[which(dat[,i]==j)]
			if(length(cur)>1){
				hashes_labels[which(labels==j)]=.hash.tip(cur, tips)
			} else {
				hashes_labels[which(labels==j)]=NA
			}
			
		}
	}
	names(hashes_labels)=labels
	hashes_labels=hashes_labels[!is.na(hashes_labels)]
	
	# redundancies for labels
	tmp=table(hashes_labels)
	if(any(tmp[tmp>1]->redund)){
		root=.hash.tip(tips, tips)
		for(r in names(redund)){
			if(r==root){
				rdx=which(hashes_labels==r)
				warning(paste("redundant labels encountered at root:\n\t", paste(names(hashes_labels[rdx]), collapse="\n\t"), sep=""))
				hashes_labels=hashes_labels[-rdx[!rdx==min(rdx)]]
			} else {
				rdx=which(hashes_labels==r)
				warning(paste("redundant labels encountered:\n\t", paste(names(hashes_labels[rdx]), collapse="\n\t"), sep=""))
				hashes_labels=hashes_labels[-rdx[!rdx==max(rdx)]]
			}
		}
	}
	
#	cat("resolving descendants for splits in tree...\n")
	tmp=hashes.phylo(phy, tips)
	hashes_tree=tmp$hash
	phy$node.label=rep("",max(phy$edge))
	mm=match(hashes_tree, hashes_labels)
	nodelabels=ifelse(is.na(mm), "", names(hashes_labels[mm]))
	nodelabels[is.na(nodelabels)]=""
	nodelabels=nodelabels[(Ntip(phy)+1):max(phy$edge)]
	tmp=table(nodelabels[nodelabels!=""])
	if(any(tmp[tmp>1]->redund)){
		for(r in names(redund)){
			rdx=which(nodelabels==r)
			nodelabels[rdx[rdx!=min(rdx)]]=""
		}
	}
	
	phy$node.label=nodelabels
	
	edges=NULL
    
	desc=.cache.descendants(phy)$tips[-c(1:Ntip(phy))]
    tidx=match(tips, phy$tip.label)
    
	FUN=function(taxon){
		nm=rownames(dat)
		dat=as.matrix(dat, ncol=ncol(dat))
		rownames(dat)=nm
		if(!taxon%in%dat) {
			try=agrep(taxon, unique(c(dat)), value=TRUE)
			if(length(try)){
				warning(paste(sQuote(taxon), " not encountered in 'taxonomy'\n\nIntended search may have been:\n\t", paste(try, collapse="\n\t", sep=""), sep=""))
			} else {
				warning(paste(sQuote(taxon), " not encountered in 'taxonomy'", sep=""))

			}
			
			return(NULL)
		}
		expected=rownames(which(dat==taxon, arr.ind=TRUE))
		if(length(expected)==1){ # single tip
			x=which(phy$tip.label==expected)
			hs=.hash.tip(expected, tips)
			res=list(unexpected=c(), missing=c())
			attr(res, "node")=x
			attr(res, "hash")=hs
			attr(res, "expected")=sort(expected)
			return(res)
			
		}
		
		bin=as.integer(tips%in%expected)
		
#		rownames(edges)=1:nrow(edges)
		N=Ntip(phy)
        
        ff=.get.parallel()

		dist=unlist(ff(desc, function(x) {
            y=tidx%in%x
            sum(abs(y-bin))
        }))
		nearest=N+which(dist==min(dist))
		hs=character(length(nearest))
		for(i in 1:length(nearest))hs[i]=.hash.tip(edges[nearest[i],], tips)
		res=lapply(nearest, function(x) {
				   tt=tips[which(tidx%in%desc[[x-N]])]
				   unexpected=sort(setdiff(tt, expected)) ## in tree but unexpected
				   missing=sort(setdiff(expected, tt)) ## expected but not in tree
				   tmp=list(unexpected=unexpected, missing=missing)
				   hs=.hash.tip(edges[x,], tips)
				   attr(tmp, "node")=x
				   attr(tmp, "hash")=hs
				   attr(tmp, "expected")=sort(expected)
				   return(tmp)
				   })
		res
	}
	
	# missed labels
	mm=match(hashes_labels, hashes_tree)
	if(any(is.na(mm))){
		mss=hashes_labels[is.na(mm)]
		if(!strict){
			N=Ntip(phy)
			
			ll=list()
			nm=rev(names(mss))
			for(x in 1:length(nm)){
				ll[[x]]=FUN(nm[x])
			}
			near_tmp=ll
			names(near_tmp)=nm
			near=lapply(1:length(near_tmp), function(idx){
						tmp=near_tmp[[idx]]
						x=nm[idx]
						if(is.null(tmp)) return(c())
						if(length(tmp)==1) {
							nd=attributes(tmp[[1]])$node
							names(nd)=paste("\"", x, "\"", sep="")
							if(phy$node.label[nd-N]=="") return(nd) else return(c())
						} else {
							return(c())
						}
			})
			near=unlist(near)
			dd=duplicated(near)
			true_missed=nm[!nm%in%gsub("\"","",names(near))]
			if(any(dd)) near=near[-which(dd)]
			phy$node.label[near-N]=names(near)
		}
	}
	
	nn=gsub("\"", "", unique(phy$node.label[phy$node.label!=""]))
	hl=names(hashes_labels)
	ms=sort(hl[!hl%in%nn])
	phy$missed=ms
	if(length(ms)){
		warning(paste("labels missing from 'phy':\n\t", paste(ms, collapse="\n\t"), sep=""))
	}

	phy$FUN=FUN

	phy
}


.TESTING.root.phylo=function(phy, outgroup, taxonomy=NULL){
## GENERAL FUNCTION FOR ROOTING (based on outgroup)
# taxonomy: classification data.frame with 'species' minimally as a rownames

	if(!is.null(sys)) {
		sys=cbind(rn=rownames(sys), taxonomy)
		rows=unique(unlist(sapply(outgroup, function(o) which(sys==o, arr.ind=TRUE)[,1])))
		outgroup=rownames(sys)[rows]
		outgroup=outgroup[outgroup%in%phy$tip.label]
	} else {
		if(!all(outgroup%in%phy$tip.label)) stop("Some 'outgroup' appear missing from 'phy'.")
	}
	
	tips=match(outgroup, phy$tip.label)
	node=getMRCA(phy,tips)
	if(node==Ntip(phy)+1){
		node=getMRCA(phy, (1:Ntip(phy))[-tips])
	}
	rooted=root(phy, node=node, resolve.root=TRUE)
	rooted
}





.TESTING.ultrametricize.phylo=function(phy, trim=c("min","max","mean","depth"), depth=NULL){
	
	phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:nrow(phy$edge)) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	
	paths=xx[1:n]
	trim=switch(match.arg(trim),
				min = min(paths),
				max = max(paths),
				mean = mean(paths),
				depth = NULL)
	
	if(is.null(trim)) {
		if(!is.null(depth)) trim=depth else stop("'depth' must be supplied if 'trim=depth'")
	}
	
	tol=diff(range(paths))
	
	cat(paste("Detected maximum difference in root-to-tip path lengths of ",tol,"\n",sep=""))
	rsc=function(phy, curdepth, depth) {phy$edge.length=phy$edge.length*(depth/curdepth); phy}

	ww=which(phy$edge[,2]<=n)
	phy$edge.length[ww]=phy$edge.length[ww]+(trim-paths[phy$edge[ww,2]])
	if(trim!=depth && !is.null(depth)) {
		phy=rsc(phy, trim, depth)
	} 
	
	if(any(phy$edge.length<0)) warning("ultrametricized 'phy' has negative branch lengths") 

	return(phy)
}


#general phylogenetic utility for returning first ancestor (as a numeric referent) of the supplied node 
#author: JM EASTMAN 2010

.get.ancestor.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,2]==node),1])
}




#general phylogenetic utility for returning all ancestors (listed as given in phy$edge[,2]) of a node 
#author: JM EASTMAN 2010

.get.ancestors.of.node <-
function(node, phy) {
	a=c()
	if(node==(Ntip(phy)+1->root)) return(NULL)
	f=.get.ancestor.of.node(node, phy)
	a=c(a,f)
	if(f>root) a=c(a, .get.ancestors.of.node(f, phy))
	return(a)
}


#general phylogenetic utility for returning the first (usually, unless a polytomy exists) two descendants of the supplied node 
#author: JM EASTMAN 2010

.get.desc.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}



#general phylogenetic utility for returning all descendants (listed as given in phy$edge[,2]) of a node (excluding the supplied node) 
#author: JM EASTMAN 2011

.get.descendants.of.node <- 
function(node, phy, tips=FALSE){
	n=Ntip(phy)
	all=ifelse(tips, FALSE, TRUE)
	out <- .Call("get_descendants", tree=list(
											  NODE = as.integer(node),
											  ROOT = as.integer(n+1),
											  ALL = as.integer(all),
											  ENDOFCLADE = as.integer(dim(phy$edge)[1]),
											  ANC = as.integer(phy$edge[,1]),
											  DES = as.integer(phy$edge[,2])),
											  PACKAGE = "geiger")
	res=out$TIPS
	if(!length(res)) res=NULL
	return(res)
} 

.mrca=function(labels, phy){
	mm=labels
	
	if(all(is.character(labels))){
		
		ll=c(phy$tip.label, phy$node.label)
		mm=match(labels, ll)
		if(any(is.na(mm))) stop("Some 'labels' not encountered in 'phy'")
		
		
	}	
	
	if(!all(is.numeric(mm))) stop("Supply 'labels' as a character or integer vector")
    
    if(length(u<-unique(mm))==1) return(u)
	
	aa=unlist(lapply(mm, function(x) .get.ancestors.of.node(x, phy)))
	tt=table(aa)
	max(as.integer(names(tt[tt==length(labels)])))
}

is.phylo=function(x) "phylo"%in%class(x)

## FUNCTIONS
## grabs most exclusive tips from clade definitions (and whose tips can then be reconciled with a taxonomic database -- a lookup table)
# allows for recursion in clade definitions (clade defined in part using another clade definition)
# returns trees representing each clade definition
# 'nested': an important variable -- this is the subset of clades that are defined (recursively) within 'clades'
phylo.clades=function(clades, phy=NULL, unplaced=TRUE){
	
## give 'phy' as a multiPhylo, named list of trees (whose labels appear in the clade defs)
## clades: 
#	clades=list(
#			 Sirenoidea=c("Siren", "Pseudobranchus"), 
#			 Ambystomatidae=c("Ambystomatidae", "Plethodontidae"), 
#			 Cryptobranchoidea=c("Hynobiidae", "Cryptobranchus_alleganiensis"),
#			 CAUDATA=c("Sirenoidea","Salamandroidea","Cryptobranchoidea")
#	)
	
    if (!is.null(phy)) {
        if (!"multiPhylo" %in% class(phy) | is.null(phynm <- names(phy))) 
            stop("Supply 'phy' as a 'multiPhylo' object with names")
        tmp = character(length(phy))
        for (i in 1:length(phy)) {
            cur = phy[[i]]
            cur$edge.length = NULL
            tre = write.tree(cur)
            tmp[i] = gsub(";", "", tre)
        }
        phy = tmp
        names(phy) = phynm
        for (i in 1:length(clades)) {
            cur = clades[[i]]
            if (any(cur == names(clades)[i])) 
                stop("Encountered self-referential clade")
            mm = match(cur, phynm)
            if (any(!is.na(mm))) {
                ww = which(!is.na(mm))
                for (j in 1:length(ww)) {
                  clades[[i]][ww[j]] = phy[[mm[ww[j]]]]
                }
            }
        }
    }
    
    ll = sapply(clades, length)
    if (any(ll == 1)) {
    	unis=clades[ll==1]
    	clades=clades[ll!=1]
        warning("Non-splitting lineages found in 'clades'")
    } else {
        unis=NULL
    }
    cc = unique(c(unlist(clades)))
    tt = cc %in% names(clades)
    nested = cc[tt]
    is.nestedclade = function(clade, nested) {
        if (clade %in% nested) 
            return(TRUE)
        else return(FALSE)
    }
    
    fetch_nestedclades = function(clade, clades, nested) {
        if (is.nestedclade(clade, nested)) {
            desc = clades[[clade]]
            new = as.list(desc)
            names(new) = desc
            tt = sapply(new, is.nestedclade, nested)
            for (i in 1:length(new)) {
                if (tt[i]) 
                  new[[i]] = fetch_nestedclades(new[[i]], clades, nested)
                else new[[i]] = NA
            }
        }
        else {
            new = NA
        }
        return(new)
    }
    
    paths_through_clades = function(clades, nested) {
        nn = names(clades)
        res = lapply(nn, function(x) {
           dd = clades[[x]]
           y = lapply(dd, fetch_nestedclades, clades, nested)
           names(y)[1:length(dd)] = dd
           y
        })
        names(res) = nn
        res
    }
    
    unplaced_phy = function(phy, cladepath) {
        if (any(names(cladepath) == "unplaced")) {
            tips = names(cladepath$unplaced)
            y = .polytomy.phylo(tips)
            y$edge.length = NULL
            new = bind.tree(phy, y)
            new = drop.tip(new, "unplaced")
            return(new)
        }
        else {
            return(phy)
        }
    }
    
    tree_cladepath = function(cladepath, nested) {
        print_group = function(cladepath) {
            xx = sapply(names(cladepath), is.nestedclade, nested)
            middle = character(length(xx))
            for (i in 1:length(xx)) {
                if (xx[i]) {
                  new = cladepath[[names(xx)[i]]]
                  middle[i] = print_group(new)
                }
                else {
                  middle[i] = names(cladepath)[i]
                }
            }
            paste("(", paste(middle, collapse = ", "), ")", sep = "")
        }
        tmp = paste(print_group(cladepath), ";", sep = "")
        phy = read.tree(text = tmp)
        return(unplaced_phy(phy, cladepath))
    }
    
    cladepaths = paths_through_clades(clades, nested)
    
    if (unplaced) {
        for (i in 1:length(cladepaths)) {
            cur = cladepaths[[i]]
            if (!is.null(unplc <- attributes(clades[[i]])$unplaced)) {
                unplaced_taxa = lapply(unplc, function(x) return(NA))
                names(unplaced_taxa) = unplc
                cladepaths[[i]]$unplaced = unplaced_taxa
            }
        }
    }
    
    alt_tip_label=function(phy, unis){
    	phy$alt.tip.label=character(Ntip(phy))
    	pt=phy$tip.label
    	for(i in 1:length(unis)){
    		if(length(cur<-unis[[i]])!=1) stop("'unis' should have a single member in each element")
    		if(!is.na(mm<-match(cur, pt))) phy$alt.tip.label[mm]=names(unis)[i]
    	}
    	phy
    }
    
    phy = lapply(cladepaths, tree_cladepath, nested)
    if(length(unis)) phy = lapply(phy, alt_tip_label, unis)
    nn=sapply(phy, Ntip)
    midx<-which(nn==max(nn))
   	master=phy[midx]
    
    ## RESOLVE NODE LABELS for 'master' tree (most comprehensive)
   	if(length(master)==1){
   		master=master[[1]]
   		tt=master$tip.label
   		null=.hash.tip(c(), tt)
		
   		mm=hashes.phylo(master, tt)
   		ss=sapply(phy, function(x) .hash.tip(x$tip.label, tt))
   		if(any(ss==null)){
			warning(paste("The following not encountered:\n\t", paste(names(ss)[which(ss==null)], collapse="\n\t")))
		}
   		ts=table(ss)
   		if(any(ts>1)){
   			drp=numeric()
   			ts=ts[ts>1]
   			for(i in 1:length(ts)){
   				ww=which(ss==names(ts[i]))
   				drp=c(drp, ww[which(ww!=min(ww))])
   			}
   			warning("The following node labels appear redundant:\n\t", paste(names(ss)[drp], collapse="\n\t"))
   			ss=ss[-drp]
   		}
   		xx=match(ss, mm$hash)
   		if(any(is.na(xx))){
   			warning(paste("The following not encountered:\n\t", paste(names(ss)[is.na(xx)], collapse="\n\t")))
   			ss=ss[-which(is.na(xx))]
   		}
   		N=Ntip(master)
   		nd=character(Nnode(master))
   		nd[xx-N]=names(ss)
   		master$node.label=nd
   		phy[[midx]]=master
   	}
    return(phy)
}


lookup.phylo=function(phy, taxonomy=NULL, clades=NULL){
## taxonomy expected to have first column at same level as tip labels in phy
## first row in taxonomy is most exclusive
## clade_defs are phylogenetic trees of a clade representation
	
	if(!is.null(taxonomy)){
		if(!any(taxonomy[,1]%in%phy$tip.label) & !is.null(rownames(taxonomy))) {
			taxonomy=as.data.frame(as.matrix(cbind(species=rownames(taxonomy), taxonomy)),stringsAsFactors=FALSE)
		}
	}
	
	related_tips=function(tips, phy){
		tips=tips[tips%in%phy$tip.label]
		if(length(tips)<2) stop("related_tips(): 'tips' to be found in 'phy' are too few")
		nd=unique(sapply(tips, function(x) match(x, phy$tip.label)))
		if(length(nd)==1) return(nd)
		anc=getMRCA(phy, nd)
		dd=.get.descendants.of.node(anc, phy, tips=TRUE)
		return(phy$tip.label[dd])
	}
	
	tips_in_group=function(phy, taxonomy=NULL, clade_def){
		
		lengths=sapply(clade_def$tip.label, function(grp){
						if(!is.null(taxonomy)){
						   ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
						   if(length(ww)>0){
						   return(TRUE)
						   } else {
						   return(FALSE)
						   } 
						} else {
							return(length(phy$tip.label[which(phy$tip.label%in%grp)])>0)
						}
		})
		
		if(sum(lengths)<2) return(NA)
		
		## FIXME: use 'lengths' to determine if 'clade_def' is satisfied (at least two members needed)
		
		
		tmp=sapply(clade_def$tip.label, function(grp){
				   if(!is.null(taxonomy)){
					ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
					if(length(ww)>0){
						return(unique(taxonomy[ww,1]))
					} else {
						return(c())
					} 
				   } else {
					return(phy$tip.label[which(phy$tip.label%in%grp)])
				   }
		})
		
		## most exclusive 'spp' (recognized from tree)
		spanning_taxa=unique(unlist(tmp))
		recovered_taxa=unique(related_tips(spanning_taxa, phy))
		
		return(recovered_taxa)
	}
	
	orig_phy=phy
	
	if(!is.null(taxonomy)){
		## check ordering of taxonomy
		oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
		if(!all(oo==c(1:ncol(taxonomy)))){
			warning("Assuming 'taxonomy' is not from most to least exclusive")
			taxonomy=taxonomy[,ncol(taxonomy):1] #reverse ordering of columns
		}
		original_taxonomy=taxonomy
		
		## check rank of tree
		phy$node.label=NULL
		gtips=sapply(phy$tip.label, function(x) unlist(strsplit(gsub(" ", "_", x), "_"))[1])
		check=sapply(list(phy$tip.label, gtips), function(x) sum(x%in%taxonomy[,1]))
		if(max(check)==check[2]) {
			genuslevel_tree=TRUE
			phy$tip.label=unname(gtips)
		} else {
			genuslevel_tree=FALSE
		}
		tips=phy$tip.label
		
		## prune taxonomy down to members in the phylogeny
		taxonomy=taxonomy[taxonomy[,1]%in%tips,]
		matching_tips=taxonomy[,1]
		
		## BUILD taxonomic mapping to each row in 'phy'
		mm=match(tips, matching_tips)
		tt=taxonomy[mm,]
		colnames(tt)=names(taxonomy)
	} else {
		tt=c()
	}
	
	if(!is.null(clades)){
		tips=phy$tip.label
		clade_defs=phylo.clades(clades)
#		cat("resolving clades...\n\t")
		res=lapply(1:length(clade_defs), function(idx) {
				   def=clade_defs[[idx]]
#				   cat(paste(names(clade_defs)[idx], "\n\t", sep=""))
				   tt=try(tips_in_group(phy, taxonomy, def),silent=TRUE)
				   if(inherits(tt, "try-error")) return(NA) else return(tt)
				   }) 
#		cat("\n")
		names(res)=names(clade_defs)
		
		zz=sapply(res, function(x) all(is.na(x)))
		if(any(zz)){
			warning(paste("taxa not represented in 'tips':\n\t", paste(names(res)[zz], collapse="\n\t"), sep=""))
		}
		res=res[!zz]
		clade_defs=clade_defs[!zz]
		
		## BUILD clade-level mapping to each row in 'phy'
		gg=matrix("", nrow=Ntip(phy), ncol=length(res))
		for(i in 1:length(res)){
			nm=names(clade_defs)[i]
			cur_tips=res[[i]]
			zz=tips%in%cur_tips
			gg[zz,i]=nm
		}
		colnames(gg)=names(clade_defs)
		tt=as.matrix(cbind(tt, gg))
	} 
	
	if(is.null(tt)){
		if(!is.null(nn<-phy$node.label)){
			names(nn)=Ntip(phy)+1:length(nn)
			nn=nn[nn!=""]
			if(length(nn)){
				tt=matrix("", nrow=Ntip(phy), ncol=length(nn))
				dd=.cache.descendants(phy)$tips
				for(i in 1:length(nn)){
					tt[phy$tip.label%in%phy$tip.label[dd[[as.integer(names(nn[i]))]]],i]=nn[i]
				}
				colnames(tt)=nn
				
			} else {
				return(NULL)
			}
		} else {
			return(NULL)
		}
	}
	phy_mapping=tt
	rownames(phy_mapping)=orig_phy$tip.label
	return(as.data.frame(as.matrix(phy_mapping), stringsAsFactors=FALSE)) 	
}



## FUNCTIONS ##

.TESTING_bind.phylo=function(phy, taxonomy){
## phy: a 'rank' level phylogeny (tips of 'phy' should be matchable to taxonomy[,rank])
## taxonomy: a mapping from genus, family, order (columns in that order); rownames are tips to be added to constraint tree
##		-- 'taxonomy' MUST absolutely be ordered from higher exclusivity to lower (e.g., genus to order)
## rank: rank at which groups are assumed to be monophyletic (currently for 'family' only)
## returns a nodelabeled constraint tree based on 'phy' and 'rank'-level constraints
	
	
	
#	oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
#	if(!all(oo==c(1:ncol(taxonomy)))){
#		warning("Assuming 'taxonomy' is not from most to least exclusive")
#		taxonomy=taxonomy[,ncol(taxonomy):1]
#	}
	taxonomy=as.data.frame(taxonomy, stringsAsFactors=FALSE)
	rank=colnames(taxonomy)[unique(which(taxonomy==phy$tip.label, arr.ind=TRUE)[,"col"])]
	if(length(rank)!=1) stop("tips in 'phy' must occur in a single column of 'taxonomy'")
	
	tax=taxonomy
	ridx=which(colnames(tax)==rank)
	tax=tax[,ridx:ncol(tax)]
	tips=rownames(tax)
	
	original_taxonomy=taxonomy
	
	
	# PRUNE 'rank'-level tree if some taxa unmatched in 'tips'
	if(any(zz<-!phy$tip.label%in%tax[,rank])){
		warning(paste("taxa not represented in 'tips':\n\t", paste(phy$tip.label[zz], collapse="\n\t"), sep=""))
		phy=drop.tip(phy, phy$tip.label[zz])
	}
	
	tmp=unique(c(phy$tip.label, phy$node.label))
	all_labels=tmp[!is.na(tmp)&tmp!=""]
	exclude=apply(tax, 1, function(x) !any(x%in%all_labels))
	missing_tips=rownames(tax)[exclude]
	if(length(missing_tips)){
		warning(paste("tips missing data in 'taxonomy' and excluded from constraint tree:\n\t", paste(missing_tips, collapse="\n\t"), sep=""))
	}
	tax=tax[!exclude,]
	
	# find tips that have data but whose data not at 'rank' level (plug in deeper in tree)
	at_rank=tax[,rank]%in%phy$tip.label
	deeper_tips=tax[!at_rank,]
	if(nrow(deeper_tips)){
		ww=apply(deeper_tips, 1, function(x) x[[min(which(x%in%all_labels))]])
		for(i in 1:length(ww)){
			nn=.nodefind.phylo(phy, ww[i])
			if(!is.null(nn)){
				tmp=.polytomy.phylo(names(ww[i]))
				phy=bind.tree(phy, tmp, where=nn)
			}
		}
		phy=compute.brlen(phy, method="Grafen")
		warning(paste("tips missing data at 'rank' level in 'taxonomy' but included in constraint tree:\n\t", paste(rownames(deeper_tips), collapse="\n\t"), sep=""))
	}
	phy$node.label=NULL
	tax=as.matrix(original_taxonomy[at_rank,1:ridx])
    if(ridx==1) {
        rownames(tax)=rownames(original_taxonomy)[at_rank]
        colnames(tax)=colnames(original_taxonomy)[ridx]
    }
    tax=as.data.frame(tax)
	if(any(is.na(tax[,rank]))) stop("Corrupted data encountered when checking taxonomy[,rank]")
	
	rsc=function(phy, age=1) {
		ee=phy$edge.length
		ag=max(heights(phy))
		phy$edge.length=ee*(age/ag)
		phy
	}
	## CREATE 'rank'-level subtrees
	mm=min(phy$edge.length)
	ss=split(tax, tax[,rank])
	subtrees=lapply(1:length(ss), function(idx) {
			  curnm=names(ss)[idx]
			  x=ss[[idx]]
			  rnm=rownames(x)
			  y=as.matrix(x)
			  d=apply(y, 2, function(z) if(all(z=="") | length(unique(z))==1) return(TRUE) else return(FALSE))
			  if(any(d)) {
				nm=colnames(y)
				y=as.matrix(y[,-which(d)])
				colnames(y)=nm[-which(d)]
			  }
			  if(!length(y)){
				cur=.polytomy.phylo(rnm, mm/2)
			  } else {
				cur=phylo.lookup(cbind(y,names(ss)[idx]))
				cur=compute.brlen(cur)
			  }
			  cur$root.edge=0
			  cur=rsc(cur, mm/2)
			  return(cur)
			  
	})
	names(subtrees)=names(ss)
	
	## PASTE in 'rank'-level subtrees
	contree=glomogram.phylo(phy, subtrees)
	contree=compute.brlen(contree)
	contree
}


.fill.taxonomy=function(taxonomy){
	push.taxon=function(taxonomy, indices, column){
		for(n in which(indices)){
			taxonomy[n,column]=paste(taxonomy[n,min((column:ncol(taxonomy))[!is.na(taxonomy[n,column:ncol(taxonomy)])])],paste("R",column,sep=""),sep=".")
		}
		return(taxonomy)
	}
	
	for(c in 1:(ncol(taxonomy)-1)){
		ii=is.na(taxonomy[,c])
		if(any(ii)) taxonomy=push.taxon(taxonomy, ii, c)
	}
	return(taxonomy)	
}

phylo.lookup=function(taxonomy) {
# GENERAL FUNCTION: convert taxonomic 'lookup' table to phylogeny
# lookup is data.frame with ranks as columns
# rowlabels of 'taxonomy' are assumed to be tips of the phylogeny
# rank corresponds to (numeric) position in rank names (R to L) to use in building taxonomic tree; if null, all ranks are used
# NOTE: taxonomic groups in 'lookup' MUST be ordered by exclusivity from R to L -- e.g., genera must precede families must precede orders
	labels=rownames(taxonomy)
    mm=nrow(taxonomy)
    tax=as.matrix(taxonomy, ncol=ncol(taxonomy))
	rownames(tax)=labels
	occurrences=table(tax)
	occurrences=occurrences[names(occurrences)!="" & occurrences>1]
    if(any(labels%in%names(occurrences))) stop(paste("The following appear to be both ranks and row.names of 'taxonomy':\n\t", paste(labels[labels%in%names(occurrences)], collapse="\n\t"), sep=""))

    if(!any(occurrences==mm)){
        tax=cbind(as.matrix(taxonomy, ncol=ncol(taxonomy)),root="root")
        rownames(tax)=labels
        occurrences=c(occurrences, root=mm)
    }
	
	env=Sys.getenv()
	if("ignoreMULTICORE"%in%names(env)) {
		f=lapply
	} else {
		if(.check.parallel()) {
			f=function(X,FUN) mclapply(X,FUN,mc.silent=TRUE)
		} else {
			f=lapply
		}
	}
	
    #clds=f(names(occurrences), function(x) {
    #		tmp=which(tax==x, arr.ind=TRUE)
	#	xcol=max(tmp[,"col"])
	#	dat=as.matrix(tax[xrow<-tmp[,"row"],1:xcol], ncol=xcol)
	#	tab=occurrences[names(occurrences)!=x]
	#	unique(sapply(1:nrow(dat), function(idx){
	#		cur=dat[idx,]
	#		   if(any(cur%in%names(tab)->xx)){
	#				ht=tab[names(tab)%in%cur]
	#				return(names(ht)[max(which(ht==max(ht)))])
	#		   } else {
	#				return(rownames(dat)[idx])
	#		   }
	#	}))
	#})
    
    
    clds=f(names(occurrences), function(x) {
		tmp=which(tax==x, arr.ind=TRUE)[,"row"]
        taxinx=rownames(tax)[tmp]
        return(taxinx)
	})
    names(clds)=names(occurrences)
    
    nn=sapply(clds, length)
    clds=clds[order(nn)]
    
    ## exclude DUPLICATES
    eq=f(1:length(clds), function(idx){
        others=clds[-which(names(clds)==names(clds)[idx])]
        vv=sapply(others, function(x) setequal(clds[[idx]], x))
        if(any(vv)) return(names(which(vv))) else return()
    })
    drop=c()
    tmp=unique(c(unlist(eq)))
    if(length(tmp)){
        mm=match(tmp, names(clds))
        tmp=tmp[order(mm)]
        mm=mm[order(mm)]
        names(mm)=tmp
        
        for(i in 1:length(mm)){
            ww=which(sapply(eq, function(x) names(mm[i])%in%x))
            tt=names(clds)[ww]
            a=which(names(clds)==names(mm[i]))
            b=match(tt, names(clds))
            if(any(b<a)) drop=c(drop, names(mm[i]))
        }
        
        clds=clds[-which(names(clds)%in%drop)]
    }
    
    resolve_definition=function(taxon, clades){
        cldother=clades[-which(names(clades)==taxon)]
        tips=clades[[taxon]]
        
        resolver=function(tips, cld){
            
            members=c()
            if(!length(cld)) {
                members=c(members, tips)
            } else {
                tt=sapply(cld, function(y){
                    all(y%in%tips)
                })
                
                if(any(tt)){
                    nm=names(cld)[max(which(tt))]
                    members=c(members, nm)
                    ftips=tips[-which(tips%in%cld[[nm]])]
                    fcld=cld[which(tt)]
                    fcld=fcld[!names(fcld)%in%nm]
                    members=c(members, resolver(ftips, fcld))
                } else {
                    members=c(members, tips)
                }

            }
            return(members)
        }
       
        def=resolver(tips, cldother)
        def
    }
    
    defs=lapply(1:length(clds), function(idx){
        taxon=names(clds)[idx]
        resolve_definition(taxon, clds)
    })
    names(defs)=names(clds)
    
    op=options()
    op$expressions=max(op$expressions, 500000)
    options(op)
	
	tmp=phylo.clades(defs)
	phy=tmp[[length(defs)]]
	tt=table(phy$tip.label)
	if(any(tt>1)){
		warning(paste("The following tips occur multiply, suggesting non-monophyly of subtending groups:\n\t", paste(names(tt[tt>1]), collapse="\n\t"), sep=""))
		warning("Offending tips removed from labeled phylogeny")
		phy=drop.tip(phy, names(tt[tt>1]))
	}
	phy$root.edge=0
	phy		
}

