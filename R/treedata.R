# Treedata is a function internal to GEIGER
# that makes sure that the names of the taxa 
# in the tree and data file match, and prunes 
# things accordingly; it returns a list with 
# two elements: phy and data

# data.names is optional, and will replace the names or rownames
# of data when matching data to the tree

# if sort is T, data will have rows in the same order
# as the taxon names in phy$tip.label

treedata<-function(phy, data, data.names=NULL, sort=FALSE, warnings=TRUE)
{

	if(is.vector(data)) {
		data<-as.matrix(data)
	}
	if(is.factor(data)) {
		data<-as.matrix(data)
	}
	if(is.array(data) & length(dim(data))==1) {
		data<-as.matrix(data)
	}
	
	if(is.null(data.names)) {
		if(is.null(rownames(data))) {
				stop("names for 'data' must be supplied")
#JME				data.names<-phy$tip.label
#JME				if(warnings)
#JME					cat("Warning: no tip labels, order assumed to be the same as in the tree\n")
		} else {
			data.names<-rownames(data)
		}
	}
	nc<-name.check(phy, data, data.names)
	if(is.na(nc[[1]][1]) | nc[[1]][1]!="OK") {
		if(length(nc[[1]]!=0)) {
			phy=drop.tip(phy, as.character(nc[[1]]))
			if(warnings) {
				warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t",
						  paste(nc[[1]], collapse="\n\t"), sep=""))
#JME			print(nc[[1]])
#JME			cat("\n")
			}
		}
	
		if(length(nc[[2]]!=0)) {
			m<-match(data.names, nc[[2]])
			data=as.matrix(data[is.na(m),])
			data.names<-data.names[is.na(m)]
			if(warnings) {
				warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t",
						  paste(nc[[2]], collapse="\n\t"), sep=""))
#JME			print(nc[[2]])
#JME			cat("\n")
			}
		}
 	}
	order<-match(data.names, phy$tip.label)	

	rownames(data)<-phy$tip.label[order]
	
	if(sort) {

    	index <- match(phy$tip.label, rownames(data))
   		data <- as.matrix(data[index,])
	}
	
	phy$node.label=NULL
	
	return(list(phy=phy, data=data))
}