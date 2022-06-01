#include <Rcpp.h>
#include <vector>

/* ANCESTORS */
/* INTERNAL C++ */
void compileancestors(int const &node,
                     int const &root,
                     int const &nrow,
                     std::vector<int> &TAXA,
                     std::vector<int> const &anc,
                     std::vector<int> const &des)
{
    // GENERAL: makes no assumption about values or ordering of 'anc' and 'des'

	/*
	 * node: current internal node
	 * nrow: number of rows in the edge matrix
	 * anc: first column of 'phylo' edge matrix (ancestors)
	 * des: second column of 'phylo' edge matrix (descendants)
	 */

	int i, d, a, eoc=nrow, r=root, n=node;

	for(i=0; i<eoc; i++) {
        d=des.at(i);
        if(d == n)
        {
            TAXA.push_back(d);
            a=anc.at(i);
            if(d != a)
            {
                if(a != r)
                {
                    compileancestors(anc.at(i), r, eoc, TAXA, anc, des);
                }
                else
                {
                    TAXA.push_back(anc.at(i));
                }
            }
		}
	}
}

/* ANCESTORS */
RcppExport SEXP compile_ancestors (SEXP mat)
{
	// using compileancestors

	try {
		/* call in parameters associated with 'phylo' object */
		Rcpp::List obj(mat);
		int node = Rcpp::as<int>(obj["node"]);
        int root = Rcpp::as<int>(obj["root"]);
        int nrow = Rcpp::as<int>(obj["nrow"]);
		std::vector<int> anc = Rcpp::as<std::vector<int> >(obj["ANC"]);
		std::vector<int> des = Rcpp::as<std::vector<int> >(obj["DES"]);

		std::vector<int> TAXA;

		/* get ancestors */
		compileancestors(node, root, nrow, TAXA, anc, des);

		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("TAXA",TAXA));

    } catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    } catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
    }
    return R_NilValue;
}


/* DESCENDANTS */
/* INTERNAL C++ */
void compiledescendants(int const &node,
                        int const &nrow,
                        std::vector<int> &TIPS,
                        std::vector<int> const &anc,
                        std::vector<int> const &des,
                        std::vector<int> const &keep)
{
    // GENERAL: makes no assumption about values or ordering of 'anc' and 'des'

	/*
	 * node: current internal node
	 * nrow: number of rows in the edge matrix
	 * TIPS: vector in which to store all descendants
	 * anc: first column of edge matrix (ancestors)
	 * des: second column of edge matrix (descendants)
     * keep: whether to continue search at element
	 */

	int i, d, eoc=nrow, n=node;

	for(i=0; i<eoc; i++) {
		if(anc.at(i) == n)
		{
            d=des.at(i);
            TIPS.push_back(d);
            if(keep.at(i)==(signed)1)
            {
                compiledescendants(d, eoc, TIPS, anc, des, keep);
            }
		}
	}
}

/* DESCENDANTS */
/* C++ | R INTERFACE; gather tip descendants main function */
RcppExport SEXP compile_descendants (SEXP mat)
{
    // using compiledescendants

	try {
		/* call in parameters associated with 'phylo' object */
		Rcpp::List phylo(mat);

		int node = Rcpp::as<int>(phylo["node"]);
		int nrow =  Rcpp::as<int>(phylo["nrow"]);
		std::vector<int> anc=phylo["ANC"];
		std::vector<int> des=phylo["DES"];
        std::vector<int> keep=phylo["keep"];

		std::vector<int> TIPS;

		/* get descendants */
		compiledescendants(node, nrow, TIPS, anc, des, keep);

		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(Rcpp::Named("TIPS",TIPS));

    } catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    } catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
    }
    return R_NilValue;
}


/* DESCENDANTS */
/* INTERNAL C++ */
void gatherdescendants(int const &node,
					   int const &root,
					   int const &endofclade,
					   std::vector<int> &TIPS,
					   std::vector<int> const &anc,
					   std::vector<int> const &des,
					   int const &all)
{
    // ASSUMES loop can terminate when reaching descendant integers (TIPS) smaller than 'root'

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

/*
 * FUNCTION TO return tip descendant given a node label and an S3 'phylo' object from the R-package ape (Paradis E, J Claude, and K Strimmer 2004 [BIOINFORMATICS 20:289-290])
 * author: Jonathan M Eastman 07.23.2011
 */

/* DESCENDANTS */
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

RcppExport SEXP cache_descendants (SEXP phy)
{
    /* requires preorder (pruningwise ordering) of 'phylo' object */
	try {
		Rcpp::List phylo(phy);
		int N = Rcpp::as<int>(phylo["N"]);
		int maxnode = Rcpp::as<int>(phylo["MAXNODE"]);
		std::vector<int> anc = Rcpp::as<std::vector<int> >(phylo["ANC"]);
		std::vector<int> des = Rcpp::as<std::vector<int> >(phylo["DES"]);

		int rows = maxnode-1;
		int root = N+1;

		std::vector< std::vector<int> > TIPS;
		std::vector< std::vector<int> > FDESC;
		std::vector< std::vector<int> > ADESC;
		std::vector< std::vector<int> > AANC;

		std::vector<int> empty;


		int i, j, k, s, t, z, dn, fd;

		/* initialize TIPS with known descendants (tips and root), otherwise leave empty */
		std::vector<int> cur;
		for(i = 0; i < maxnode; i++) {
			FDESC.push_back(empty);
            AANC.push_back(empty);
			if(i < N)
			{
				cur.push_back(i+1);
				TIPS.push_back(cur);
				ADESC.push_back(cur);
			}
			else
			{
				TIPS.push_back(empty);
				ADESC.push_back(empty);
			}
			cur.clear();
		}

		/* store nodes associated with root -- TIPS */
		for(i=0; i<N; i++){
			cur.push_back(i+1);
		}
		TIPS.at(N)=cur;

		/* store nodes associated with root -- ALL */
		cur.clear();
		for(i=0; i<maxnode; i++){
			if(i!=N)
			{
				cur.push_back(i+1);
			}
		}
		ADESC.at(N)=cur;

		/* store nodes associated with root -- FIRST */
		cur.clear();
		for(i=0; i<rows; i++){
			int idx = anc.at(i);
			if(idx==root)
			{
				cur.push_back(des.at(i));
			}
		}
		FDESC.at(N)=cur;

		/* collect descendants of each node in edge matrix (using pruningwise order to eliminate unnecessary computations) */
		for(i = 0; i < rows; i++){
			int nd = des.at(i);
			if(nd > N) {
				/* find descendant nodes */
				std::vector<int> subtends;
				for(j=0; j<rows; j++){
					int idx = anc.at(j);
					if(idx==nd)
					{
						subtends.push_back(des.at(j));
					}
				}
				FDESC.at(nd-1)=subtends;

				/* find immediate descendants of nd */
				std::vector<int> subtendedtips;
				std::vector<int> subtendednodes;
				s=subtends.size();
				for(k = 0; k < s; k++){

					/* find nodes subtended by immediate descendants of nd */
					fd = subtends.at(k);
					subtendednodes.push_back(fd);

					if(fd<root) {
						subtendedtips.push_back(fd);
					} else {
						std::vector<int> descnodes = ADESC[fd-1];
						t=descnodes.size();
						for(z = 0; z < t; z++){
							dn=descnodes[z];
							if(dn<root)
							{
								subtendedtips.push_back(dn);
							}
							subtendednodes.push_back(dn);
						}
					}
				}
				/* store tips associated with nd into main list */
				TIPS.at(nd-1)=subtendedtips;
				ADESC.at(nd-1)=subtendednodes;
                s=subtendednodes.size();

                /* store ancestors */
                for(k=0; k<s; k++){
                    int idx = subtendednodes.at(k);
                    std::vector<int> ancnodes = AANC[idx-1];
                    ancnodes.push_back(nd);
                    AANC.at(idx-1)=ancnodes;
                }
			}
		}

		for(i=0; i<N; i++){
			ADESC.at(i)=empty;
		}
        for(k=0; k<maxnode; k++){
            int idx = k+1;
            if(idx!=root){
                std::vector<int> ancnodes = AANC[k];
                ancnodes.push_back(root);
                AANC.at(k)=ancnodes;
            }
        }
		return Rcpp::List::create(Rcpp::Named("tips",TIPS),
								  Rcpp::Named("fdesc",FDESC),
								  Rcpp::Named("adesc",ADESC),
                                  Rcpp::Named("anc",AANC));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
	}
	return R_NilValue;
}
