`area.between.curves` <-
function(x, f1, f2, xrange=c(0,1))
{
	a<-0.0;
	for(i in 1:length(x)) {
		if(x[i]>=xrange[1] & x[i]<=xrange[2]) {
			if(i==1) {
				lhs<-0  
			} else if(x[i-1]<xrange[1]) {
					lhs<-xrange[1]
				} else lhs<-x[i-1];
			if(i==length(x)) { 
				rhs<-x[i]
			} else if(x[i+1]>xrange[2]) {
				rhs<-xrange[2]; 
			} else rhs<-x[i+1];
			a<-a+(f2[i]-f1[i])*(rhs-lhs)/2;	
		} else if(i!=1) if(x[i-1]>=xrange[1] & x[i-1]<=xrange[2]) {
			y1<-f1[i-1]+(f1[i]-f1[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
			y2<-f2[i-1]+(f2[i]-f2[i-1])*(xrange[2]-x[i-1])/(x[i]-x[i-1])
			a<-a+(y2-y1)*(xrange[2]-x[i-1])/2;
		} else if(i!=length(x)) if(x[i+1]>=xrange[1] & x[i+1]<=xrange[2]) {
			y1<-f1[i]+(f1[i+1]-f1[i])*(xrange[1]-x[i])/(x[i+1]-x[i])
			y2<-f2[i]+(f2[i+1]-f2[i])*(xrange[1]-x[i])/(x[i+1]-x[i])

			a<-a+(y2-y1)*(x[i+1]-xrange[1])/2;
		}
	}
	return(a)
}

