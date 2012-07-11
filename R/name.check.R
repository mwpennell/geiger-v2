`name.check` <-
function(phy, data)
{
	
#	if(is.null(data.names)) 
#	{
	if(is.vector(data)){
		data.names<-names(data)
	} else {
		data.names<-rownames(data)
	}
#	}
	t<-phy$tip.label
	r1<-t[is.na(match(t,data.names))]
	r2<-data.names[is.na(match(data.names,t))]
	
	r<-list(sort(r1), sort(r2))
	
	names(r)<-cbind("tree_not_data", "data_not_tree")
	if(length(r1)==0 && length(r2)==0) return("OK")
	else return(r)
}

