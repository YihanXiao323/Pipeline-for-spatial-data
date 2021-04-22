countmatrix=$1
locationfile=$2
resultfilename=$3
DEmethod=$4


# check input param
if [ $# -lt 4 ]; then
  echo "usage:\n Integrate.sh <countmatrix> <locationfile> <resultfilename> <Demehtod (Trendsceek, Binspect)>"
  exit
fi


if [ ! -d ${resultfilename} ]; then
  mkdir ${resultfilename}
fi

#Preprcessing
Rscript Preprocessing/Giotto_preprocess.r -r ${resultfilename} -c ${countmatrix} -l ${locationfile} 

#DE genes
if [ $DEmethod = "Binspect" ]; then
      Rscript DEgenes/Binspect.R -r ${resultfilename}
elif [ $DEmethod = "Trendsceek" ]; then
      Rscript DEgenes/Trendsceek.R -r ${resultfilename}
fi

#Clustering
Rscript Clustering/HMRF.r -r ${resultfilename}
#Rscript Clustering/Kmeansclustering.R -r ${resultfilename}
