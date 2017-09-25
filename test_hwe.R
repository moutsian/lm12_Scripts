#function to calculate p-value for deviations from HWE using a chi square test
calc_hwe=function(hom1,het,hom2){

n=hom1+het+hom2
pA=(2*hom1+het)/(2*n)
pa=1-pA
observed=c(hom1,het,hom2)
expected=n*c(pA^2,2*pA*pa,pa^2)
res=sum(((observed-expected)^2)/expected)
pval=1-pchisq(res,df=1)
return(pval)
}
#END