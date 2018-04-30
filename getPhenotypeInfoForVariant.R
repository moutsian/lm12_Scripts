getPheno=function(chr,pos,build=38) {
    if( build == 38) {
        server="http://rest.ensembl.org"
    } else if( build == 37) {
        server="http://grch37.rest.ensembl.org"
    } else {
        print("ERROR: Invalid build supplied.")
    }

    ext=paste0("/phenotype/region/homo_sapiens/",chr,":",pos,"-",pos,"?feature_type=Variation")
    r=GET(paste(server, ext, sep = ""), content_type("application/json"))
    stop_for_status(r)
    p=fromJSON(toJSON(content(r), null="null"))
    if(length(p)==0) {
        ext=paste0("/phenotype/region/homo_sapiens/",chr,":",pos-1e5,"-",pos+1e5,"?feature_type=Variation")
        r=GET(paste(server, ext, sep = ""), content_type("application/json"))
        stop_for_status(r)
        p=fromJSON(toJSON(content(r), null="null"))
        l=p$phenotype_associations
        names(l)=p$id
        l2=mapply(function(x,y) {cbind(x,"id"=y)}, l, names(l))
        pheno=rbindlist(lapply(l2,flatten),fill=T)
        return(pheno)
    } else {
        l=p$phenotype_associations
        names(l)=p$id
        l2=mapply(function(x,y) {cbind(x,"id"=y)}, l, names(l))
        pheno=rbindlist(lapply(l2,flatten),fill=T)
        return(pheno)
    }

}


library(httr,lib="/software/team152/Rpackages/") 
library(jsonlite,lib="/software/team152/Rpackages/")
library(data.table,lib="/software/team152/Rpackages/")
library(data.table,lib="/software/team152/Rpackages/")
 library(curl,lib="/software/team152/Rpackages/")
getPheno(2,60585806)
