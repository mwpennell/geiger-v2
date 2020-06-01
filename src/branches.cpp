#include <Rcpp.h>
#include <cmath>
#include <vector>

/* C++ | R INTERFACE; pruning algorithm of BM likelihood (based on diversitree:::make.bm.direct */
RcppExport SEXP bm_direct (SEXP dat, SEXP pars)
{
	/*
	 * all objects ordered from 1:(Nnode(phy)+Ntip(phy)) unless noted otherwise
	 * dat: a list of elements
		* len: edge lengths (vector)
		* root: root ID
		* y: tip data
		* order: order of internal IDs for pruning algorithm
	 * pars: rates associated with each branch
	 */

	try {
		/* call in parameters associated with 'dat' object */
		Rcpp::List cache(dat);

		int root = Rcpp::as<int>(cache["root"]);
		int n = Rcpp::as<int>(cache["n"]);

        double drift = Rcpp::as<double>(cache["drift"]);

		std::vector<double> len = Rcpp::as<std::vector<double> >(cache["len"]);
		std::vector<double> y = Rcpp::as<std::vector<double> >(cache["y"]);
		std::vector<double> var = Rcpp::as<std::vector<double> >(cache["var"]);
        std::vector<double> known = Rcpp::as<std::vector<double> >(cache["given"]); //
		std::vector<int> intorder = Rcpp::as<std::vector<int> >(cache["intorder"]);
		std::vector<int> tiporder = Rcpp::as<std::vector<int> >(cache["tiporder"]);
		std::vector<int> descR = Rcpp::as<std::vector<int> >(cache["descRight"]);
		std::vector<int> descL = Rcpp::as<std::vector<int> >(cache["descLeft"]);

		std::vector<double> rates = Rcpp::as<std::vector<double> > (pars);

		std::vector<double> lq;
		lq.assign(n,0.0);

		double yi, ri, li, m1, m2, v1, v2, v12, m12, n12, nm, nv, m, mm, v, k;

		double const PIx = 4.0*atan(1.0);

		std::vector<double> branchinitM;
		std::vector<double> branchinitV;
		std::vector<double> branchbaseM;
		std::vector<double> branchbaseV;

		branchinitM.assign(n,0.0);
		branchinitV.assign(n,0.0);
		branchbaseM.assign(n,0.0);
		branchbaseV.assign(n,0.0);

		int i, z, cur, d1, d2;

		/* mean and variance for leaves */
		z=tiporder.size();
		for(i=0; i<z; i++){
			cur=tiporder[i]-1;
			yi=y[cur];
			li=len[cur];
			ri=rates[cur];

			branchinitM[cur] = yi;
			branchbaseM[cur] = yi + drift*li;
			branchbaseV[cur] = var[cur] + li*ri;
		}

		/* mean, variance, and density for edges */
		z=intorder.size();
		for(i=0; i<z; i++){
			cur=intorder[i]-1;
			d1=descR[cur]-1;
			d2=descL[cur]-1;
            m1=branchbaseM[d1];
			m2=branchbaseM[d2];
			v1=branchbaseV[d1];
			v2=branchbaseV[d2];

            v12=v1+v2;

            m = (((m1*v2) + (m2*v1))/v12);                  // phylogenetic mean expectation
            branchinitM[cur] = m;

            v = ((v1*v2)/v12);
            branchinitV[cur] = v;

            m12=pow((m1-m2),2);
            lq[cur] = ((-m12/(2*v12)) - (log(2*PIx*v12)/2));

            k=known[cur];

            if( k == (signed)1 )                            // resolve whether node state is given (k==1)
            {
                nm=y[cur];
                nv=var[cur];
                mm=m;

                v12=v+nv;
                m = ((mm*nv) + (nm*v))/v12;
                branchinitM[cur] = m;

                v = (v*nv)/v12;
                branchinitV[cur] = v;

                m12=pow((mm-nm),2);
                lq[cur]+=((-m12/(2*v12)) - (log(2*PIx*v12)/2));
            }

            li=len[cur];
            branchbaseM[cur] = m + drift*li;
            branchbaseV[cur] = v + rates[cur]*li;
		}

		/* compute root */
		cur=root-1;
		d1=descR[cur]-1;
		d2=descL[cur]-1;
		m1=branchbaseM[d1];
		m2=branchbaseM[d2];
		v1=branchbaseV[d1];
		v2=branchbaseV[d2];
		v12=v1+v2;

        m=(((m1*v2) + (m2*v1))/v12); //
		branchinitM[cur] = m; //
        v=((v1*v2)/v12); //
		branchinitV[cur] = v;
		m12=pow((m1-m2),2);
		lq[cur] = ((-m12/(2*v12)) - (log(2*PIx*v12)/2));

        /* compute root lnL (either ML or given) */
        k=known[cur];
        if(k == (signed)1 )
        { // given state
            nm=y[cur];
            nv=var[cur] + v;
            n12=pow((nm-m),2);
            lq[cur] += ((-n12/(2*nv)) - (log(2*PIx*nv)/2));
        }
        else
        { // ML state
            lq[cur] += (- (log(2*PIx*v)/2));

        }


		/* PREPARE OUTPUT FOR R */
		return Rcpp::List::create(
								  Rcpp::Named("initM",branchinitM),
								  Rcpp::Named("initV",branchinitV),
								  Rcpp::Named("baseM",branchbaseM),
								  Rcpp::Named("baseV",branchbaseV),
								  Rcpp::Named("lq",lq)
								  );


    } catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    } catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
    }
    return R_NilValue;
}


/* C++ | R INTERFACE; determine which branches subtended by node are not in list of excluded nodes (and their descendants)  */
RcppExport SEXP open_subtree (SEXP dat, SEXP desc)
{
	/*
	 * dat: a list of elements
	 *  node: of interest
	 *  exclude: vector of excluded nodes
	 *  N: tips
	 * desc: a list of ALL descendants from 1:(Ntip(phy)+Nnode(phy))
	 */

	try {
		/* call in parameters associated with 'dat' object */
		Rcpp::List cache(dat);
		Rcpp::List adesc(desc);
		int N = Rcpp::as<int>(cache["N"]);
		int n = Rcpp::as<int>(cache["n"]);
		int nn = N+n;
		int node = Rcpp::as<int>(cache["node"]);
		std::vector<int> exclude = Rcpp::as<std::vector<int> >(cache["exclude"]);
		std::vector<int> subtended;
		std::vector<int> drop;

		int i, j, k, se, sn, sd, nd, cur;

		std::vector<int> res;

		if(node>N)
		{
			/* fix exclude vector (drop) */
			se=exclude.size();
			for(i=0; i<se; i++){
				cur=exclude[i];
				drop.push_back(cur);
				if( (cur>node) & (cur<=nn) )
				{
					std::vector<int> dd=adesc[cur-1];
					sd=dd.size();
					for(j=0; j<sd; j++){
						drop.push_back(dd[j]);
					}
				}
			}

			std::vector<int> nodedesc = adesc[node-1];
			se=drop.size();
			if(se>0)
			{
				/* find descendants not in drop list */
				sn=nodedesc.size();
				sd=drop.size();
				for(i=0; i<sn; i++){
					cur=nodedesc[i];
					k=0;
					for(j=0; j<sd; j++){
						nd=drop[j];
						if(cur!=nd)
						{
							k+=1;
						}
					}
					if(k==sd)
					{
						res.push_back(cur);
					}
				}
			}
			else
			{
				sn=nodedesc.size();
				for(i=0; i<sn; i++){
					res.push_back(nodedesc[i]);
				}
			}
		}

		/* PREPARE OUTPUT FOR R */
		return Rcpp::wrap(res);

    } catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    } catch(...) {
		::Rf_error( "C++ exception: unknown reason" );
    }
    return R_NilValue;
}
