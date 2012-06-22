#include <Rcpp.h>
#include <vector>
/* 
 * GENERAL FUNCTIONS for dealing with an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 01.11.2011
 */

/* INTERNAL C++ */
void sortedges(std::vector<double> const &unsortededges,
			   std::vector<double> &edges,
			   std::vector<int> const &des)
{
	std::vector<int>::size_type i,j;
	//int i, j;
	for(i=0; i<edges.size(); i++) {
		for(j=0; j<des.size(); j++){
			if(des.at(j)==(signed)i+1)
			{
				edges.at(i)=unsortededges.at(j);							 
			}
		}
	}
}

/* INTERNAL C++ */
void gatherdescendants(int const &node,
					   int const &root,
					   int const &endofclade,
					   std::vector<int> &TIPS,
					   std::vector<int> const &anc, 
					   std::vector<int> const &des,
					   int const &all)
{
	/* 
	 * node: current internal node
	 * root: most internal node
	 * endofclade: number of rows in the phylogenetic edge matrix (of the 'phylo' object)
	 * TIPS: vector in which to store all (tip) descendants of node (by node label)
	 * anc: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
	 * des: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
	 * all: whether to return all (all=1) or just tips (all=0)
	 */
	
	int i, eoc=endofclade, r=root, n=node, aa=all;

	for(i=0; i<eoc; i++) {
		if(anc.at(i) == n)
		{
			if(des.at(i)>r)
			{
				if(all==1)
				{
					TIPS.push_back(des.at(i));
				}
				
				gatherdescendants(des.at(i), r, eoc, TIPS, anc, des, aa);
				
			}
			else 
			{
				TIPS.push_back(des.at(i));
	
			}
		}
	}
}
	
					   
/* INTERNAL C++ */
void descend_vcv(int const &node,
				 double const &edge,
				 int const &root,
				 int const &endofclade,
				 std::vector<int> const &anc, 
				 std::vector<int> const &des, 
				 std::vector<double> &vcv)
{
	/* 
	 * node: internal node of which descendants are recorded
	 * edge: length to add to covariances of all descendants where 'node' is internal
	 * root: node to initiate vcv update
	 * endofclade: rows in the phylogenetic edge matrix (of the 'phylo' object)
	 * anc: first column of phylo edge matrix (e.g., phy$edge[,1]
	 * des: second column of phylo edge matrix (e.g., phy$edge[,2]
	 * vcv: the current variance-covariance matrix
	 */
	
	int n;
	std::vector<int>::size_type a,b,j;
	n=root-1;
	std::vector<int> TIPS;	
	TIPS.reserve(root-1);
	
	gatherdescendants(node, root, endofclade, TIPS, anc, des, 0);
		
	/* add 'edge' length of 'node' to covariances V[a,b] where a and b are tips descended from 'node' and a!=b [OFFDIAGONAL ELEMENTS in V] */
	for(a=0; a<TIPS.size(); a++) {
		for(b=a; b<TIPS.size(); b++) {
			if(a!=b)
			{
				vcv[ TIPS.at(a)-1 + n*( TIPS.at(b)-1 ) ] = vcv[ TIPS.at(b)-1 + n*(TIPS.at(a)-1 ) ] += edge;
			}
		}
	}
	
	/* add 'edge' length of 'node' to variances of V[a,b] where a and b are tips descended from 'node' and a==b [DIAGONAL ELEMENTS in V]*/ 
	for(j=0; j<TIPS.size(); j++) {
		vcv[ TIPS.at(j)-1 + n*( TIPS.at(j)-1 ) ] += edge;
	}
	
	std::vector<int>().swap(TIPS);
}


/* INTERNAL C++ */
void vcv_internal (int const &maxnode,
				   int const &root,
				   int const &endofclade,
				   std::vector<int> const &anc, 
				   std::vector<int> const &des, 
				   std::vector<double> const &edges, 
				   std::vector<double> &vcv)
{
	/* 
	 * maxnode: largest internal node
	 * root: most internal node
	 * endofclade: number of rows in the phylogenetic edge matrix (of the 'phylo' object)
	 * anc: first column of phylo edge matrix (e.g., phy$edge[,1]
	 * des: second column of phylo edge matrix (e.g., phy$edge[,2]
	 * edges: sorted branch lengths (by node label), including the root (ntips+1)
	 * vcv: the current state of the variance-covariance matrix 
	 */
	
	int n, nn, r=root, mx=maxnode, eoc=endofclade;	
	nn=r-1;
	
	/* cycle over internal nodes within the tree; tree stem does not contribute to VCV so is skipped */
	for(n=r+1; n<=mx; n++) {
		
		/* find descendants of n; add edge length of n to variances for tips and covariances among tips */	
		descend_vcv(n, edges.at(n-1), r, eoc, anc, des, vcv);
		
	}
	
	/* add leaves to lengths in V[n,n] where n is a tip [DIAGONAL ELEMENTS in V]*/
	for(n=1; n<r; n++) {
		vcv[ n-1 + nn*( n-1 ) ] += edges.at(n-1);
	}
}



/* 
 * FUNCTION TO GENERATE VARIANCE-COVARIANCE MATRIX from an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 01.11.2011
 */

/* C++ | R INTERFACE; variance-covariance main function */
RcppExport SEXP vmat (SEXP tree) 
{
	/* 
	 * tree: a list of elements 
		* ROOT: most internal node
		* MAXNODE: least internal internal node (and largest valued in edge matrix)
		* ENDOFCLADE: rows in edge matrix
		* ANC: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
		* DES: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
		* EDGES: edge lengths, sorted by node label and including the root (at position Ntip(phy)+1)
		* VCV: the current state of the variance-covariance matrix, initialized with 0s in R
	 */
		
	try {
		std::vector<int>::size_type i;
		
		/* call in parameters associated with 'phylo' object */
		Rcpp::List phylo(tree);

		int root = (int) Rcpp::as<int>(phylo["ROOT"]);
		int maxnode =  Rcpp::as<int>(phylo["MAXNODE"]);
		int endofclade =  Rcpp::as<int>(phylo["ENDOFCLADE"]);
		std::vector<int> anc=phylo["ANC"];
		std::vector<int> des=phylo["DES"];
		std::vector<double> unsortededges=phylo["EDGES"];
		std::vector<double> edges=phylo["EDGES"];
		std::vector<double> V=phylo["VCV"];
		
		/* initialize edges */
		for(i=0; i<edges.size(); i++) {
			edges.at(i)=0;
		}		
		
		/* sort edges by node label */
		sortedges(unsortededges, edges, des);
		
		/* call to primary function that updates VCV matrix */
		vcv_internal(maxnode, root, endofclade, anc, des, edges, V);
		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("VCV",V));


    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}


/* 
 * FUNCTION TO return tip descendant given a node label and an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 07.23.2011
 */

/* C++ | R INTERFACE; gather tip descendants main function */
RcppExport SEXP get_descendants (SEXP tree) 
{
	/* 
	 * tree: a list of elements 
	 * NODE: node for which descendants will be returned
	 * ROOT: most internal node
	 * ALL: whether to gather all (ALL=1) or just tips (ALL=0)
	 * ENDOFCLADE: rows in edge matrix
	 * ANC: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
	 * DES: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
	 */
	
	try {		
		/* call in parameters associated with 'phylo' object */
		Rcpp::List phylo(tree);
		
		int node = Rcpp::as<int>(phylo["NODE"]);
		int root = Rcpp::as<int>(phylo["ROOT"]);
		int all = Rcpp::as<int>(phylo["ALL"]);
		int endofclade =  Rcpp::as<int>(phylo["ENDOFCLADE"]);
		std::vector<int> anc=phylo["ANC"];
		std::vector<int> des=phylo["DES"];
		
		std::vector<int> TIPS;	
		if(all==0)
		{
			TIPS.reserve(root-1);
		}
		else 
		{
			TIPS.reserve(2*(root-1));	
		}
		
		/* get descendants */
		gatherdescendants(node, root, endofclade, TIPS, anc, des, all);

		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("TIPS",TIPS));
	
    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}

/* 
 * FUNCTION TO GENERATE binary representation of all edges in tree from an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 07.24.2011
 */

/* C++ | R INTERFACE; main function */
RcppExport SEXP binary_edges (SEXP tree) 
{
	/* 
	 * tree: a list of elements 
	 * NTIP: number of species in tree
	 * ROOT: most internal node
	 * ENDOFCLADE: rows in edge matrix
	 * MAXNODE: least internal internal node
	 * ANC: first column of 'phylo' edge matrix (e.g., phy$edge[,1])
	 * DES: second column of 'phylo' edge matrix (e.g., phy$edge[,2])
	 * EDGES: initialized binary matrix
	 */
	
	try {		
		/* call in parameters associated with 'phylo' object */
		Rcpp::List phylo(tree);
		
		int ntip = Rcpp::as<int>(phylo["NTIP"]);
		int root = Rcpp::as<int>(phylo["ROOT"]);
		int endofclade =  Rcpp::as<int>(phylo["ENDOFCLADE"]);
		int maxnode = Rcpp::as<int>(phylo["MAXNODE"]);
		std::vector<int> anc=phylo["ANC"];
		std::vector<int> des=phylo["DES"];
		std::vector<int> edges=phylo["EDGES"];

		int i,j;
		std::vector<int>::size_type k;
		int cur_node, cur_tip, cur_spp;
		
		std::vector<int> TIPS;	
		
		for(i=0; i<maxnode; i++)
			for(j=0; j<ntip; j++)
				edges[(ntip*i)+j]=0;
			
		
		/* collect tip descendants of each node */
		for(i=0; i<maxnode; i++){
			cur_node=i+1;
			if(cur_node<root)
			{
				cur_tip=cur_node;
				for(j=0; j<ntip; j++){
					cur_spp=j+1;
					if(cur_tip==cur_spp)
					{
						edges[(ntip*i)+j]=1;
					}
				}
			}
			else 
			{
				TIPS.reserve(root-1);
				gatherdescendants(cur_node, root, endofclade, TIPS, anc, des, 0);
				for(k=0; k<TIPS.size(); k++) {
					cur_tip=TIPS.at(k);
					for(j=0; j<ntip; j++){
						cur_spp=j+1;
						if(cur_tip==cur_spp)
						{
							edges[(ntip*i)+j]=1;
						}
					}
				}
				std::vector<int>().swap(TIPS);
			}
		}
		
		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("EDGES",edges));
		
    } catch( std::exception &ex ) {		
		forward_exception_to_r( ex );
    } catch(...) { 
		::Rf_error( "C++ exception: unknown reason" ); 
    }
    return R_NilValue; 
}


