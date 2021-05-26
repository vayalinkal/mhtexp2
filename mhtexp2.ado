program mhtexp2
    version 14
    syntax varlist [if] [in], treatment(varlist) [ subgroup(varname) combo(string) exclude(name) only(name) bootstrap(integer 3000) controls(varlist) treatnames(string) studentized(integer 1) idbootmat(name)]

    if ("`combo'" != "" & "`combo'" != "pairwise" & "`combo'" != "treatmentcontrol"){
        display "INVALID combo choose either pairwise or treatmentcontrol"
        error
    }

	global ncontrols: word count `controls'  
	
    if ("`exclude'" == "") mata: excludemat = (.,.,.)
    else mata: excludemat = `exclude'
    if ("`only'" == "") mata: onlymat = (.,.,.)
    else mata: onlymat = `only'

    mata: Y = buildY("`varlist'")
    mata: D = buildD("`treatment'")
    mata: DX = buildDX("`treatment' `controls'")

	if ($ncontrols > 0){ 
		mata: X = DX[.,2::cols(DX)]
	}
	else {
		mata: X = 0
	}
    mata: sub = buildsub("`subgroup'", D)
    mata: sizes = buildsizes(Y, D, sub)
    mata: combo = buildcombo("`combo'", sizes[3])
	mata: outrows = tokens("`treatnames'")
    mata: numpc = buildnumpc(combo)
    mata: select = buildselect(onlymat, excludemat, sizes[1], sizes[2], numpc)
		if ("`idbootmat'" == ""){
		mata: rseed(0)
		mata: idbootmat = runiformint(rows(Y), `bootstrap', 1, rows(Y)) // an n by B matrix of simulated samples of all the units with replacement
	} 
    mata: results = seidelxu2(Y, sub, D, combo, select, `bootstrap', DX, X, `studentized', idbootmat)
    mata: buildoutput("results", results, outrows)

    matlist results
end

mata:

    function buildY(string scalar outcomes){
        Y = st_data(., tokens(outcomes))
        return(Y)
    }
    function buildD(string scalar treatment){
        D = st_data(., tokens(treatment))
        return(D)
    }
    function buildX(string scalar controls){
        X = st_data(., tokens(controls))
        return(X)
    }
    function buildDX(string scalar indvars){
        DX = st_data(., tokens(indvars))
        return(DX)
    }
    function buildsub(string scalar subgroup, real matrix D){
        if (subgroup == ""){
            sub = J(rows(D), 1,1)
        }else{
            sub = st_data(., (subgroup))
        }
        return(sub)
    }
    function buildsizes(real matrix Y, real matrix D, real matrix sub){
        numoc = cols(Y)
        numsub = colnonmissing(uniqrows(sub))
        numg = rows(uniqrows(D)) - 1

        return((numoc, numsub, numg))
    }
    function buildcombo(string scalar strcombo, real scalar numg){
        if (strcombo == "pairwise"){
    		combo = nchoosek((0::numg), 2)
    	}else{
    		combo = (J(numg,1,0), (1::numg))
    	}
        return(combo)
    }
    function buildnumpc(real matrix combo){
        return(rows(combo))
    }
    function buildselect(real matrix only, real matrix exclude, real scalar numoc, real scalar numsub, real scalar numpc){
        if (rownonmissing(only) != 0){
            select = mdarray((numoc, numsub, numpc),0)
            for (r = 1; r <= rows(only); r++){
                i = only[r, 1]
                j = only[r, 2]
                k = only[r, 3]
                put(1, select, (i,j,k))
            }
        }else{
            select = mdarray((numoc, numsub, numpc), 1)
        }
        if (rownonmissing(exclude) !=0){
            for (r=1; r <= rows(exclude); r++){
                i = exclude[r, 1]
                j = exclude[r, 2]
                k = exclude[r, 3]
                put(0, select, (i,j,k))
            }
        }
        return(select)
    }

end


mata:

function seidelxu2(Y, sub, D, combo, select, bootstrap, DX, X, studentized, idbootmat){


// parameters set by the function
n = rows(Y)  // the number of units
// n
B = bootstrap // number of simulated samples
numoc = cols(Y) // number of outcomes
numsub = colnonmissing(uniqrows(sub)) // number of subgroups
numg=rows(uniqrows(D)) - 1 // the number of treatment groups (not including the control group)
numpc=rows(combo) // the number of pairs of treament and control groups of interest
stud = studentized // turn studentization on/off
// pi_prob = treatprobs // treatment probabilities for studentization


Nact= mdarray((numoc, numsub, numg+1), 0) // sample size of actual data for all H
regact = mdarray((numoc, numsub, numpc), .) // coefficient of treatment
abregact = mdarray((numoc, numsub, numpc), .) // absolute value of coefficient of treatment


for (i=1; i <= numoc; i++)
{
	for (j=1; j<=numsub; j++)
	{ 
		if(cols(DX) > 1) {
			sg = (sub :== j)
			cursgX = select(X, sg)
			barXz = mean(cursgX)
			varXz = quadvariance(cursgX)
		}
		else { 
			sg = (sub :== j)
			barXz = 0
			varXz = 0
		}

		for (l=1; l <= numpc; l++)
		{
			w = (sub :== j :& (D :== combo[l,1] :| D :== combo[l,2]))
			pi_2 = sum(sub :== j :& D :== combo[l,2])/n //sum(w)
			pi_1 = sum(sub :== j :& D :== combo[l,1])/n //sum(w)
			pi_z = sum(sub :== j)/n //sum(w)
// 			pi_t = pi_prob[l,2]
// 			pi_s = pi_prob[l,1]
			curDX = select(DX, w) // subdata of interest
			curD = select(D, w) // subdata of interest
			curY = select(Y[.,i],w) // outcomes of interest

			curDX[.,1] = (curD :== combo[l,2]) // treatment dummy

			regres = runreg(curDX, curY, barXz, varXz, pi_2, pi_1, pi_z, colsum(sg))  

			put(regres[1], regact, (i,j,l))
			put(abs(regres[1])/(stud*regres[2] + (1-stud)), abregact, (i,j,l))
// 			put(abs(regres[1])/(stud*regres[3] + (1-stud)), sturegact, (i,j,l))
// 			printf("og")
// 			abs(regres[1])/(stud*regres[2] + (1-stud))
// 			regres[2]

		}
	}
}


idboot =  idbootmat
// idboot
statsboot= mdarray((B, numoc, numsub, numpc), 0) // test statistics for all the simulated samples
regboot= mdarray((B, numoc, numsub, numpc), 0) // beta1s for all the simulated samples
abregboot= mdarray((B, numoc, numsub, numpc), 0) // |beta1|s for all the simulated samples
meanboot = mdarray((numoc, numsub, numg+1), 0) // samples means of the simulated sample for all H
varboot = mdarray((numoc, numsub, numg+1), 0) // sample variance of the simulated sample for all H
Nboot = mdarray((numoc, numsub, numg+1), 0) // sample sizes of simulated sample for all H
diffboot = mdarray((numoc, numsub, numpc),.) // difference in means for each treatment control pairs

sprintf("bootstrap iteration")
for (i=1; i <= B; i++)
{
	if (mod(i,10) == 0) printf("-")
	if (mod(i,1000) == 0) printf("%g\n", i)
	Yboot = Y[idboot[.,i], .] // all outcomes for ith simulated sample
	subboot = sub[idboot[.,i], .] // all subgroup ids for the ith simulated sample
	Dboot = D[idboot[.,i], .] // all treatment control status for the ith simulated sample
	DXboot = DX[idboot[.,i], .] // all data for the ith simulated sample
	if(cols(DXboot) > 1) {
		Xboot = X[idboot[.,i], .] // covariate data for the ith simulated sample
	}
	mdregarr = mdarray((numoc, numsub, rows(combo)), .)
	mdabregarr = mdarray((numoc, numsub, rows(combo)), .)
	mdsturegarr = mdarray((numoc, numsub, rows(combo)), .)
	for (j=1; j <= numoc; j++)
	{
		for (k=1; k <= numsub; k++)
		{		
			if(cols(DXboot) > 1) {
				sg = (subboot :== k)
				cursgX = select(Xboot, sg)
				barXz = mean(cursgX)
				varXz = quadvariance(cursgX)
			}
			else { 
				sg = (sub :== k)
				barXz = 0
				varXz = 0
			}
			for (l=0; l <= numg; l++)
			{
				w = (subboot :== k :& Dboot :== l)
				put(mean(Yboot[.,j], w), meanboot, (j, k, l+1))
				put(variance(Yboot[.,j], w), varboot, (j,k,l+1))
				CP = quadcross(w, 0, Yboot[.,j] , 1)
				put(CP[cols(CP)], Nboot, (j, k, l+1))
			}
			for (l=1; l <= numpc; l++)
			{
				w = (subboot :== k :& (Dboot :== combo[l,1] :| Dboot :== combo[l,2]))
				pi_2 = sum(sub :== k :& Dboot :== combo[l,2])/n //sum(w)
				pi_1 = sum(sub :== k :& Dboot :== combo[l,1])/n //sum(w)
				pi_z = sum(sub :== k)/n
// 				pi_t = pi_prob[l,2]
// 				pi_s = pi_prob[l,1]
				curDX = select(DXboot, w) // subdata of interest
				curD = select(Dboot, w) // subdata of interest
				curY = select(Yboot[.,j],w) // outcomes of interest
				curDX[.,1] = (curD :== combo[l,2]) // convert the treatment column of curDX to a dummy
// 				printf("boot %f", i)
// 				abs(regres[1])
// 				abs(regres[1]-get(regact, (j,k,l)))
// 				get(regact, (j,k,l))
				regres = runreg(curDX, curY, barXz, varXz, pi_2, pi_1, pi_z, colsum(sg)) // run it
				put(abs(regres[1]-get(regact, (j,k,l)))/(stud*regres[2]+(1-stud)), abregboot, (i,j,k,l))
			}
		}
	}

}

pact = mdarray((numoc,numsub, numpc), 0 ) // a matrix of 1-p values of the actual data
pboot = mdarray((B, numoc, numsub, numpc), 0) // a matrix of 1-p values of the simulated data

for (i=1; i<=numoc; i++)
{
	for (j=1; j<=numsub; j++)
	{
		for (k=1; k<=numpc; k++)
		{
// 			get(abregboot,(.,i,j,k))
			p = 1 - (sum(get(abregboot,(.,i,j,k)) :>= get(abregact, (i,j,k)) * J(B,1,1))) / B
// 			printf("pval:%f",p)
			put(p, pact, (i,j,k))
			for (l=1; l<=B; l++)
			{
				sp = 1 - (sum(get(abregboot, (.,i,j,k)) :>= get(abregboot, (l,i,j,k)) * J(B,1,1))) / B;
				put(sp, pboot, (l,i,j,k))
			}
		}
	}
}

/* calculate p-values based on single hypothesis testing */
alphasin = mdarray((numoc, numsub, numpc), 0) // the smalled alpha's that reject H based on single testing procedure (Remark 3.2)

for (i=1; i<=numoc; i++)
{
	for (j=1; j<=numsub; j++)
	{
		for (k=1; k<=numpc; k++)
		{
			ptemp =  get(pboot, (.,i,j,k))
			sortp = sort(ptemp, -1)
			v = (get(pact, (i,j,k)) * J(B,1,1 )) :>= sortp
			indx = find(v)
			if (indx == NULL) {
			    q = 1
			}else{
			    q = indx/B
			}
			put(q, alphasin, (i,j,k))
		}
	}
}

psin = mdarray((numoc, numsub, numpc), 0) // p-values based on single hypothesis testing (psin = alphasin)
for (k=1; k<=numpc; k++){
	put(get(alphasin, (.,.,k)), psin, (.,.,k))
}

/* Calculate p-values based on multiple hypothesis tesitng */
nh = 0 // the number of hypothesis
for (k=1; k <= numpc; k++)
{
	nh = nh + sum((*(select[k,.]))[.,.])
}

statsall = J(nh, 8+B, 0)
// columns 1-5 present the id's of the hypotheses, outcomes, subgroups, and treatment (control) groups;
// the 6th column shows the studentized differences in means for all the hypotheses based on the actual data
// the 7th column presents p-values based on single hypothesis testing;
// the 8th column presents 1-p values based on the actual data;
// the subsequent columns present the corresponding 1-p values based on the simulated samples

counter=1
for (i=1; i<=numoc; i++)
{
	for (j=1; j<=numsub; j++)
	{
		for (k=1; k<=numpc; k++)
		{
			if ( (*(select[k, .]))[i,j] == 1 ){
				rowvect = (i,j,k)
				// abs removed here
				statsall[counter, .] = (counter, i, j, combo[k, .], get(regact, rowvect) , get(psin, rowvect), get(pact, rowvect), get(pboot, (., i, j, k))')
				counter = counter + 1
			}
		}
	}
}
statsrank = sort(statsall, 7) // rank the rows according to the p-values based on single hypothesis testin
alphamul = J(nh, 1, 0) // the smallest alpaha that reject the hypothesis based on Theorem 3.1
alphamulm = J(nh, 1, 0) // the smallest alpha's that reject the hypothesis based on Remark 3.8


for (i=1; i<=nh; i++)
{
	maxstats = colmax(statsrank[(i::rows(statsrank)), (9::cols(statsrank))]) //maximums of 1-p values for all remaining H for all simulated samples
	sortmaxstats = sort(maxstats', -1)'
	v = statsrank[i, 8] :>= sortmaxstats
	indx = find(v)
	if (indx == NULL){
	    q = 1
	}else{
	    q = indx/B
	}
	alphamul[i] = q
	if (i==1){
		alphamulm[i]=alphamul[i]
	}else{
		sortmaxstatsm=J(1,B,0) // compute at each quantile the maximum of critical values for all the "true" subst of H
		for (j=nh-i+1; j >= 1; j--)
		{
			subset = nchoosek(statsrank[(i::rows(statsrank)), 1], j) // all the subsets of H with j elements
			sumcont = 0 // total number of subsets of H with j elements that contradict any of previously rejected H
			for (k=1; k<=rows(subset); k++ ){
				cont = 0 // cont = 1 if any of the previously rejected hypothesis contradicts the current subset of H
				for (l=1; l <= i-1; l++)
				{
					tempA = statsall[(subset[k,.]), (2..3)]
					tempB = J(rows(tempA), 1, statsrank[l, (2..3)] )
					sameocsub = select(subset[k,.], (ismember(tempA, tempB,1)')) // the H with same outcome as the lth H
					if (cols(sameocsub) >= 1){
						tran = mat2cell(statsall[(sameocsub), (4..5)], J(1, cols(sameocsub),1) , 2) // cell array that presents sets of equal treatment(control) groups implied by "transitivity" under the null H in sameocsub
						trantemp=tran
					}
					if (cols(sameocsub) <= 1){
						cont = 0
						maxstatsm = colmax(statsall[(subset[k,.]), (9::cols(statsall))]) // maximums of the 1-p values within the subset of H for all the simulated samples
						sortmaxstatsm = colmax(sortmaxstatsm \ sort(maxstatsm', -1)')
						break
					}else{
						counter = 1
						while ( max(asarray_keys(tran)[.,1]) > max(asarray_keys(trantemp)[.,1]) || counter == 1 ){
							tran=trantemp
							trantemp = asarray_create("real", 2)
							asarray(trantemp, (1,1), asarray(tran, (1,1)))
							counter=counter + 1
							for (m=2; m<=max(asarray_keys(tran)[.,1]); m++){
								belong = 0 // total number of rows "transtemp" that "tranm" can be connected to by "transivity"
								for (n=1; n <= max(asarray_keys(trantemp)[.,1]); n++){
									trantempn = asarray(trantemp, (n,1))
									tranm = asarray(tran, (m,1))
									unq = uniqrows( (trantempn, tranm)' )'
									test = unq :< cols(trantempn) + cols(tranm)
									if (sum(test) == cols(unq)){
										asarray(trantemp, (n,1), unq)
										belong = belong+1
										if (n==max(asarray_keys(trantemp)[.,1]) && belong ==0){
											asarray(trantemp, (n+1,1), tranm)
										}
									}
								}
							}
						}
						for (p=1; p<=max(asarray_keys(tran)[.,1]); p++){
							if (sum(ismember(statsrank[l, (4..5)], asarray(tran, (p,1)), 2)) == 2){ // the lth previously rejected H contract the current subset of H
								cont=1
								break
							}
						}
					}
					if (cont==1){
						break
					}
				}
				sumcont=sumcont+cont
				if (cont==0){
					maxstatsm = colmax(statsall[(subset[k,.]), (9::cols(statsall))])
					sortmaxstatsm = colmax(sortmaxstatsm \ sort(maxstatsm', -1)')
				}
			}
			if (sumcont==0){
				break; // if all the subsets of H with j elements do not contradict any of the previously rejected hypothesis, smaller subsets do not either
			}
		}
		indx = find(statsrank[i,8] :>= sortmaxstatsm)
		if (indx == NULL){
			qm = 1
		}else{
			qm = indx/B
		}
		alphamulm[i] = qm
	}
}

bon = rowmin((statsrank[.,7]*nh, J(nh,1,1) )) // p-values based on the Bonferroni method
holm = rowmin((statsrank[.,7]:*(nh::1), J(nh,1,1))) // p-values based on the Holm's method

output = sort((statsrank[.,(1::7)], alphamul, alphamulm, bon, holm),1)
output = output[., (2::cols(output))]

return(output)

}

end

// functions 

mata:

function runreg(matrix X, colvector y, rowvector Xbar, matrix Xvar, scalar pi_2, scalar pi_1, scalar pi_z, scalar full_n)
{
	n    = rows(X)
	D = X[., 1] // treatment column
	Ddummy = (D :== 1)
	Cdummy = (D :== 0)
	y1 = select(y, Ddummy)
	y0 = select(y, Cdummy)
	if(cols(X) > 1){
		DX = X
		X = X[.,2::cols(X)]  // covariates
		X1 = select(X, Ddummy)
		X0 = select(X, Cdummy)
		DX1 = (J(rows(select(D, Ddummy)),1,1), X1)
		DX0 = (J(rows(select(D, Cdummy)),1,1), X0)
	} 
	else {
		DX = D
		X1 = (J(rows(select(D, Ddummy)),1,0))
		X0 = (J(rows(select(D, Cdummy)),1,0))
		DX1 = (J(rows(select(D, Ddummy)),1,1))
		DX0 = (J(rows(select(D, Cdummy)),1,1))
	}
	
	// run the regressions
	DX1pDX1  = quadcross(DX1, DX1)
	DX1pDX1i = invsym(DX1pDX1)
	b1 = DX1pDX1i*quadcross(DX1, y1) 

	DX0pDX0  = quadcross(DX0, DX0)
	DX0pDX0i = invsym(DX0pDX0)
	b0 = DX0pDX0i*quadcross(DX0, y0)
	
	if(cols(X) > 1){
		bX1 = b1[2::rows(b1),.]
		e1 = y1 - X1*bX1
		s1 = quadvariance(e1) 
		
		bX0 = b0[2::rows(b0),.]
		e0 = y0 - X0*bX0
		s0 = quadvariance(e0)
	}
	else {
		s1 = quadvariance(y1) 
		s0 = quadvariance(y0)
		bX1 = 0
		bX0 = 0
	}
	
	
	
	est_ATE = (b1[1,1]-b0[1,1]) + Xbar*(bX1-bX0)
// 	b1
// 	b0
// 	(b1-b0)
// 	(est_ATE)
// 	((b1-b0)+ Xbar*(bX1-bX0))
// 	((b1-b0)\(bX1-bX0))
	
	est_VAR = (1/pi_2)*(s1) + (1/(pi_1))*s0 + (1/(pi_z))*((bX1 - bX0)'*Xvar*(bX1 - bX0))
// 	est_SE = sqrt((est_VAR*full_n)/(full_n-cols(DX)))
	est_SE = sqrt(est_VAR)
//  	est_SE = 1
// est_SE = sqrt((est_VAR*full_n))
		
	return (est_ATE, est_SE)
}

function mdarray(transmorphic rowvector rowvec , fill)
{
	if (cols(rowvec) == 3) c_l = 1
	else c_l = rowvec[1,4]
	r_k = rowvec[1,3]


	a = J(r_k,c_l,NULL)
	for (k=1; k<=rows(a); k++)
	{
		for (l=1; l<=cols(a); l++)
		{
			a[k,l] = &J(rowvec[1,1],rowvec[1,2],fill)
		}
	}
	return(a)
}

function put(val,matrix x, rowvector rowvec)
{
	/* Usage: value to put, matrix to put it in, i,j of dimension k, to put it at.*/

	if (cols(rowvec)== 3) c_l = 1
	else c_l = rowvec[1,4]
	r_k = rowvec[1,3]
	i = rowvec[1,1]
	j = rowvec[1,2]

	(*(x[r_k,c_l]))[i,j]=val
}

function get(matrix x, rowvector rowvec)
{
	/* Usage: matrix to get from, i,j of dimension k, of value to get. */


	if (cols(rowvec)== 3) c_l = 1
	else c_l = rowvec[1,4]
	r_k = rowvec[1,3]
	i = rowvec[1,1]
	j = rowvec[1,2]

	return((*(x[r_k,c_l]))[i,j])
}

function nchoosek(colvector V, real scalar K)
{
    A = J(comb(rows(V), K), K, .)
    com = J(100, 1, .)
    n = rows(V)
    for (i = 1; i <= K; i++)
    {
        com[i] = i
    }
    indx = 1
    while (com[K] <= n ){
        for (i = 1; i <= K; i++)
        {
            A[indx,i] = V[com[i]]
        }
        indx = indx+1
   

        t = K
        while (t != 1 && com[t] == n - K + t)
        {
            t = t - 1
        }
        com[t] = com[t] + 1;
        for (i = t +1; i <= K; i++)
        {
            com[i] = com[i-1] + 1
        }
    }

    return(A)
}

transmorphic scalar function find(transmorphic vector V)
{
    // FInds the first nonzero index of a col vector V
    indx = NULL
	if (rows(V) <= 1){
		counter = cols(V)
	}else{
		counter = rows(V)
	}
    for (i=1; i <= counter; i++){
	if (V[i] != 0){
	    indx = i
	    break
	}
    }

    return(indx)
}

real matrix function ismember(real matrix A, real matrix B, real scalar r){
    //For array and A and B of same number of cols
    //returns an array of the same size as
    // A where A is in B = 1 0 otherwise
    // r == 1 compares rows.  If r == 0, then A and B should be rowvectors and it will return
	// it will return (1 x cols(A)) or cols(B)) whichever is smaller
    // where result[i,j] = 1 if A[i,j] == B[i,j] else 0
    // NOTE: r==1 is only supported when A and B are same size
    if (r == 1){
        res = J(rows(A), 1, .)
        for (i = 1; i <= rows(A); i++){
            if (sum(A[i,.] :== B[i,.]) == cols(A)) res[i] = 1
            else res[i] = 0
        }
    }else{
		if (cols(A) < cols(B)){
			small = A
			large = B
		}else{
			small = B
			large = A
		}
		res = J(cols(small), 1, .)
		for (i = 1; i <= cols(small); i++){
			test = small[i] :== large
			if (sum(test) >= 1 ){
				res[i] = 1
			}else{
				res[i] = 0
			}
		}
    }
    return(res)
}

function mat2cell(transmorphic matrix A, transmorphic vector rowD, transmorphic vector colD){
	/// the sum of rowD must equal the cols(A)
	/// the sum of colD must equal rows(A)
	rowcount = 1
	colcount = 1
	matcell = asarray_create("real", 2)
	for (i = 1; i <= cols(rowD); i++){
		for (j = 1; j <= cols(colD); j++){
			cell = A[rowcount::rowcount + rowD[i]-1, colcount..colcount + colD[j] -1]
			asarray(matcell, (i,j), cell)
			colcount = cols(cell) + 1
		}
		colcount = 1
		rowcount = rowcount + rows(cell)
	}

	return(matcell)
}

void function buildoutput(string scalar name, real matrix output, string matrix rows){
    headers = ("outcome","subgroup","treatment1","treatment2","coefficient","Remark3_2","Thm3_1", "Remark3_8", "Bonf","Holm")
    blanks = J(cols(headers), 1, "")

    headersmatrix = (blanks, headers')

    st_matrix(name, output)
    st_matrixcolstripe(name, headersmatrix)
	
	if(cols(rows) == rows(output)) {
    blanks2 = J(rows(output), 1, "")
    rowsmatrix = (blanks2, rows')
	
    st_matrixrowstripe(name, rowsmatrix)
	}
	else if(cols(rows) != rows(output) & cols(rows) > 0) {
	printf("treatnames length does not match number of treatments!")
	}
}

end

