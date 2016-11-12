"************************************************************************
*		Projet : Classification non supervisée			                        *
*-----------------------------------------------------------------------*
*	Kévin HOARAU		|	M1 Informatique			                                *
*	32001076	    	|	Univ. Réunion	                                      *
*************************************************************************"



"******************************
*  Déclaration des fonctions  *
*******************************"


#-------------------------------------------------------------------------------#
# 		K-Means Direct                         				                            #
#-------------------------------------------------------------------------------#
# Entrée : - data, le jeu de données                                       	    #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
#	   - up_centroid, fonction de calcul des centroides			                      #
#	   - k, nombre de classe que l'on souhaite obtenir                    	      #
#	   - centroids, matrice des centroide initiaux (aléatoire si non renseigné)   #
#	   - max_iter, nombre d'iteration maximum (illimité si non renseigné)    			#
# Sorties : - result$clusters, liste des classes obtenues                       #
#	   - result$silhouette, silhouette de la partiton                     	      #
#	   - result$inertie, inertie intra classe de la partition                     #
#	   - result$iter, nombre d'iteration                                   		    #
#-------------------------------------------------------------------------------#
kmeans_direct <- function(data,disimilarity,up_centroid,k,centroids=NULL,max_iter=0){
  nb_elem = dim(data)[1] # Nombre d'objets
  nb_var = dim(data)[2] # Nombre de variables
  centroid_change=TRUE
  iter = 0
  
  if(nb_elem>=k){
    if(is.null(centroids)){
      centroids = data[sample(nb_elem,k),]	# On choisit k centroides au hasard
    }
    
    # Debut du K-means
    while(centroid_change && (iter < max_iter || max_iter==0)){
      centroid_change=FALSE 
      mat_disim = disimilarity(data,centroids) # Calcul de la matrice des disimilarite
      clusters = rep(list(matrix(nrow=0,ncol=nb_var)), k) # On initialise une partition vide de k classes

      # Pour chaque objets
      for(i in 1:nb_elem){
      	cluster = which.min(mat_disim[i,]) # Determine le centroide le plus proche
      	clusters[[cluster]] = rbind(clusters[[cluster]],data[i,]) # Ajoute l'objet à la classe correspondante
      }
      
      # Pour chaque centroides
      for(i in 1:length(clusters)){
      	old_centroid = centroids[i,]
      	centroids[i,] = up_centroid(clusters[[i]]) # Mise à jour du centroide
      	
      	if(!identical(old_centroid,centroids[i,])){ # Detecte si le centroide à changé
      	  centroid_change=TRUE
      	}
      }
      iter = iter + 1
    }
    # Fin du K-means
  }
  else{
    clusters = list(data)
  }
  
  result <- list()
  result$clusters = clusters
  result$silhouette = silhouette(clusters,disimilarity)
  result$inertie = inertie_intra(clusters,disimilarity,up_centroid)
  result$iter = iter
  
  return(result)
}


#-------------------------------------------------------------------------------#
# 		K-Means Incremental                         				                      #
#-------------------------------------------------------------------------------#
# Entrée : - data, le jeu de données                                       	    #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
#	   - up_centroid, fonction de calcul des centroides			                      #
#	   - k, nombre de classe que l'on souhaite obtenir                    	      #
#	   - centroids, matrice des centroide initiaux (aléatoire si non renseigné)   #
#	   - max_iter, nombre d'iteration maximum (illimité si non renseigné)    			#
# Sorties : - result$clusters, liste des classes obtenues                       #
#	   - result$silhouette, silhouette de la partiton                     	      #
#	   - result$inertie, inertie intra classe de la partition                     #
#	   - result$iter, nombre d'iteration                                   		    #
#-------------------------------------------------------------------------------#
kmeans_incremental <- function(data,disimilarity,up_centroid,k,centroids=NULL,max_iter=0){

  nb_elem = dim(data)[1] # Nombre d'objets
  nb_var = dim(data)[2] # Nombre de variables
  clusters = rep(list(matrix(nrow=0,ncol=nb_var)), k) # On initialise une partition vide de k classes
  centroid_change=TRUE
  cluster_affectation = matrix(rep(0,nb_elem*2),nrow=nb_elem,ncol=2) # Matrice d'affection des objets
  iter = 0
  if(is.null(centroids)){
    centroids = data[sample(nb_elem,k),] # On choisit k centroides au hasard
   }
  
  # Debut du K-means
  while(centroid_change && (iter < max_iter || max_iter==0)){
    centroid_change=FALSE 
    
    # Pour chaque objets
    for(i in 1:nb_elem){
    
      mat_disim = disimilarity(data,centroids) # Calcul de la matrice des disimilarite
      cluster = which.min(mat_disim[i,]) # Determine le centroide le plus proche
      old_cluster = cluster_affectation[i,1] # Ancienne classe de l'objet
      pos_in_cluster = cluster_affectation[i,2] # Position dans cette classe
      
      if(old_cluster!=cluster){ # Si la classe de l'objet à changé
     
      	if(old_cluster!=0){ # Si l'objet se trouvais dans une autre classe
      	  clusters[[old_cluster]] = clusters[[old_cluster]][-pos_in_cluster,] # Enleve l'objet de la classe
      	  old_centroid = centroids[old_centroid,]
      	  centroids[old_centroid,] = up_centroid(clusters[[old_centroid]]) # Mise à jour du centroide de cette classe
      	  if(!identical(old_centroid,centroids[old_centroid,])){ # Detecte si le centroide à changé
      	    centroid_change=TRUE
      	  }
      	}
      	
      	clusters[[cluster]] = rbind(clusters[[cluster]],data[i,]) # Ajoute l'objet à sa nouvelle classe correspondante
      	cluster_affectation[i,1] = cluster # Declare la classe dans laquelle se trouve l'objet
      	cluster_affectation[i,2] = dim(clusters[[cluster]])[1] # Declare la position de l'objet dans cette classe
      	
      	old_centroid = centroids[cluster,]
      	centroids[cluster,] = up_centroid(clusters[[cluster]]) # Mise à jour du centroide de la nouvelle classe
      	if(!identical(old_centroid,centroids[cluster,])){ # Detecte si le centroide à changé
      	  centroid_change=TRUE
      	}
      }
    }
    iter = iter + 1
  }
  # Fin du K-means
  
  result <- list()
  result$clusters = clusters
  result$silhouette = silhouette(clusters,disimilarity)
  result$inertie = inertie_intra(clusters,disimilarity,up_centroid)
  result$iter = iter
  
  return(result)
}


#-------------------------------------------------------------------------------#
# 		Bisecting K-Means                         				                        #
#-------------------------------------------------------------------------------#
# Entrée : - data, le jeu de données                                       	    #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
#	   - up_centroid, fonction de calcul des centroides			                      #
#	   - sse, fonction de calcul du SSE de 2 classes                              #
#	   - k, nombre de classe que l'on souhaite obtenir                    	      #
#	   - nb_iter, nombre d'iteration du bisecting  (10 si non renseigné)    			#
# Sorties : - result$clusters, liste des classes obtenues                       #
#	   - result$silhouette, silhouette de la partiton                     	      #
#	   - result$inertie, inertie intra classe de la partition                     #
#	   - result$iter, nombre d'iteration                                   		    #
#-------------------------------------------------------------------------------#
bisecting_kmeans <- function(data,disimilarity,up_centroid,sse,k,nb_iter=10){
  clusters = list(data) # Intialise la liste des classes avec une classe contenant tout les objets
  iter=0
  nb_cluster = length(clusters) # nombre de classe
  
  while(k>nb_cluster){ # tant que le nombre de classe en est inferieur a k
    selected = iter%%nb_cluster+1  # On selectionne une classe 
    cluster = clusters[[selected]] 
    lowest_sse = -1 
    
    for(i in 1:nb_iter){ # On execute le nombre d'iteration souhaité
      sub_clusters = kmeans_direct(cluster,disimilarity,up_centroid,2)$clusters # Découpe la classe en 2 sous classe
      sub_clusters_sse = sse(sub_clusters,disimilarity,up_centroid) # calcul le SSE de ces classes
      if(lowest_sse>sub_clusters_sse || lowest_sse<0){ # Stocke ces classes si elle on un SSE plus faible que les précedantes
      	lowest_sse_sub_clusters = sub_clusters
      	lowest_sse = sub_clusters_sse
      }
    }
    
    clusters[[selected]] = lowest_sse_sub_clusters[[1]] # Remplace dans la liste la classe par la premier sous classe
    if(length(lowest_sse_sub_clusters)>1){
      clusters = c(clusters,list(lowest_sse_sub_clusters[[2]])) # Ajoute la deuxieme sous classe a la liste
    }
    nb_cluster = length(clusters)
    iter = iter +1
  }
  
  result <- list()
  result$clusters = clusters
  result$silhouette = silhouette(clusters,disimilarity)
  result$inertie = inertie_intra(clusters,disimilarity,up_centroid)
  result$iter = iter
  
  return(result)
}


#-------------------------------------------------------------------------------#
# 		Fonction de contruction du Minimum Spaning Tree			                      #
#-------------------------------------------------------------------------------#
# Entrée : - data, le jeu de données                                       	    #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
# Sorties : - le Minimum Spaning Tree                                           #
#-------------------------------------------------------------------------------#
build_mst <- function(data,disimilarity){
  nb_elem = dim(data)[1] # Nombre d'objets
  nb_var = dim(data)[2]
  free_points = t(rbind(seq(nb_elem),t(data))) # Crée une matrice indexé des données, contenant la liste des points qui ne sont pas dans l'arbre
  start_point = free_points[sample(nb_elem,1),] # Prend un point au hasard pour la racine de l'arbre
  free_points = free_points[-start_point[1],] # Enleve le point racine des données
  
  mst = matrix(ncol=3,nrow=0) # Initialise la matrice contenant l'arbre
  points = matrix(start_point,ncol=nb_var+1) # Place le points racine dans la matrices contenant la liste de point deja dans l'arbre
  
  
  while(length(free_points)>0){ # Tant que des points ne sont pas dans l'arbre
    
    free_points = matrix(data=free_points,ncol=nb_var+1) # Formate la matrice de donnée indexé
    disim = disimilarity(matrix(free_points[,-1],ncol=2),matrix(points[,-1],ncol=2)) # Disimilarite entre les points de l'arbre et les points n'y étant pas
    
    closest = which.min(disim)-1 # Intersection des points les proche proche tel que le premier est dans l'arbre et le second non 
    closest_col = closest%/%dim(free_points)[1]+1
    closest_row = closest%%dim(free_points)[1]+1
    
    q = points[closest_col,][1] # Index du point dans l'arbre
    p = free_points[closest_row,][1] # Index du point qui n'est pas dans l'arbre
    d = disim[closest_row,closest_col] # Disimilarité entre ces points
    mst = rbind(mst,c(q,p,d)) # Ajoute ces données a l'arbre
    
    points = rbind(points,free_points[closest_row,]) # Ajoute le second point à la liste des points deja dans l'arbre
    free_points = free_points[-closest_row,] # Enleve le second point de la liste des points qui ne sont pas dans l'arbre
  }
  return(mst)
}


#-------------------------------------------------------------------------------#
# 		Fonction de changement de classe d'un sous arbre d'un MST                 #
#-------------------------------------------------------------------------------#
# Entrée : - index, position de la racine du sous arbre                    	    #
#	   - clusters, liste de classes	                                              #
#	   - from_c, position de la classe de depart des elements dans clusters	      #
#	   - from_c, position de la classe d'arrivé des elements dans clusters	      #
#	   - mst, minimum spaning tree du jeu de données                      	      #
# Sorties : - clusters, liste de classes avec les elements déplacés             #
#-------------------------------------------------------------------------------#
change_cluster_mst <- function(index,clusters,from_c,to_c,mst){
  pos_in_cluster = which(clusters[[from_c]][,1]==index) # Positon de l'element a deplace dans l'ancienne classe
  clusters[[to_c]] = rbind(clusters[[to_c]],clusters[[from_c]][pos_in_cluster,]) # Ajoute l'element a la nouvelle classe
  clusters[[from_c]]= clusters[[from_c]][-pos_in_cluster,] # Enleve l'element de l'ancienne classe

  childrens = mst[,2][mst[,1]==index] # Liste des enfants de l'element
  for(i in childrens){ # Pour chaque enfants
   clusters = change_cluster_mst(i,clusters,from_c,to_c,mst) # Deplace l'enfant dans la nouvelle classe
  }
  return(clusters)
}


#-------------------------------------------------------------------------------#
# 		MST-based divisive                          				                      #
#-------------------------------------------------------------------------------#
# Entrée : - data, le jeu de données                                       	    #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
#	   - up_centroid, fonction de calcul des centroides			                      #
#	   - k, nombre de classe que l'on souhaite obtenir                    	      #
# Sorties : - result$clusters, liste des classes obtenues                       #
#	   - result$silhouette, silhouette de la partiton                     	      #
#	   - result$inertie, inertie intra classe de la partition                     #
#-------------------------------------------------------------------------------#
mst_based_divisive <- function(data,disimilarity,up_centroid,k){
  nb_elem = dim(data)[1] # Nombre d'objets
  nb_var = dim(data)[2]
  mst = build_mst(data,disimilarity) # Construit le MST
  clusters = list(t(rbind(seq(nb_elem),t(data)))) # Intialise la liste des classes avec une classe contenant tout les objets et ajoute des indexes
  cluster_affectation = rep(1,nb_elem) # Matrice d'affection des objets
  
  nb_cluster = length(clusters) # nombre de classe
  
  while(k>nb_cluster){ # Tant que l'on a moin de k classe

    pos_in_mst = which.max(mst[,3]) # Position dans la matrice du MST de la plus grande disimilarité
    index = mst[pos_in_mst,][2] # Element a partir du quel on doit crée la nouvelle classe
    mst = mst[-pos_in_mst,] # Coupe la branche
    mst = matrix(mst,ncol=nb_var+1)
    
    new_cluster = matrix(ncol=nb_var+1,nrow=0) # Initialise la nouvelle classe
    clusters = c(clusters,list(new_cluster)) # L'ajoute a la lite des classe
    
    from_cluster = cluster_affectation[index] # Position de l'ancienne classe des éléments dans la liste des classes
    to_cluster = length(clusters) # Position de la nouvelle classe des éléments dans la liste des classes
    
    clusters=change_cluster_mst(index,clusters,from_cluster,to_cluster,mst) # Change de classe l'element et ses enfants
    for(i in 1:dim(clusters[[to_cluster]])[1]){ # Met a jour les affectations des objets
      cluster_affectation[clusters[[to_cluster]][i,1]]=to_cluster
    }
    
    nb_cluster = length(clusters) # nombre de classe
  }
  
  for(i in 1:nb_cluster){ # Enleve les indexes dans les classes
   clusters[[i]] = matrix(clusters[[i]],ncol=3)[,-1]
  }
  
  result <- list()
  result$clusters = clusters
  result$silhouette = silhouette(clusters,disimilarity)
  result$inertie = inertie_intra(clusters,disimilarity,up_centroid)
  
  return(result)
}


#-------------------------------------------------------------------------------#
# 		Fuzzy K-Means                                				                      #
#-------------------------------------------------------------------------------#
# Entrée : - data, le jeu de données                                       	    #
#	   - m, indice de flou			                                                  #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
#	   - termination, seuil de tolerance pour l'arret du fuzzy K-Means            #
#	   - k, nombre de classe que l'on souhaite obtenir                    	      #
#	   - maxit, nombre d'iteration maximum (10 si non renseigné)    		        	#
# Sorties : - result$fuzzy_clusters, partition floue obtenue                    #
#	   - result$iter, nombre d'iteration                                   		    #
#-------------------------------------------------------------------------------#
fuzzy_kmeans <- function(data,m,disimilarity,termination,k,maxit=10){
  nb_elem = dim(data)[1] # Nombre d'objets
  nb_var = dim(data)[2] # Nombre de variables
  iter = 0
  centroids = matrix(nrow=k,ncol=nb_var) # Initialisation de la matrice des centroides
  
  #Partition initiale
  fuzzy_clusters = matrix(data=rep(0,k*nb_elem),nrow=k,ncol=nb_elem) 
  for (i in 1:nb_elem) {
    for (j in 1:(k-1)) {
      fuzzy_clusters[j,i]=sample(seq(0,1-sum(fuzzy_clusters[,i]),0.1),1) 
    }
    fuzzy_clusters[k,i]=1-sum(fuzzy_clusters[,i])
  }
  difference = termination +1
  
  #Debut Fuzzy K-Means
  while(difference > termination && iter < maxit){
    fuzzy_clusters_old = fuzzy_clusters
    
    #Calcul des centroides
    for (i in 1:k) {
      for (j in 1:nb_var) {
        z = fuzzy_clusters[i,]^m
        centroids[i,j] = sum(z*data[,j])/sum(z)
      }
    }
    
    # Mise a jour de la partition
    disim = disimilarity(data,centroids)
    for (i in 1:nb_elem) {
      if(sum(disim[i,]==0)<1){ # Si aucun element a une disimilarité nulle avec un centroide
        for (j in 1:k) {
          power = 2/(m-1)
          fuzzy_clusters[j,i] = 1/sum((disim[i,j]/disim[i,])^power)
        }
      }
      else{ # Sinon
        zeros = which(disim[i,]==0)
        l_zeros = length(zeros)
        fuzzy_clusters[,i] = rep(0,k)
        if(l_zeros>1){
          for (j in 1:(l_zeros-1)) {
            fuzzy_clusters[zeros[j],i]=sample(seq(0,1-sum(fuzzy_clusters[,i]),0.1),1)
          }
        }
        fuzzy_clusters[zeros[l_zeros],i]=1-sum(fuzzy_clusters[,i])
      }
    }
    difference = norm(fuzzy_clusters-fuzzy_clusters_old)
    iter = iter + 1
  }
  #Fin Fuzzy K-Means
  
  result <- list()
  result$fuzzy_clusters = fuzzy_clusters
  result$iter = iter
  
  return(result)
}


#-------------------------------------------------------------------------------#
# 		Fontion de calcul du coefficient de silhouette d'une partition            #
#-------------------------------------------------------------------------------#
# Entrée : - clusters, partition dont on souhaite calculé la silhouette    	    #
#	         - disim, fonction de generation de la matrice de disimilarite	      #
# Sortie : - coefficient de silhouette de la partition	                        #
#-------------------------------------------------------------------------------#
silhouette <- function(clusters,disim){
  iter = 0
  s = 0
  for (i in 1:length(clusters)) {
    
    for (j in 1:(dim(clusters[[i]])[1])) {
      point = matrix(clusters[[i]][j,],ncol=2)
      a = 0
      b = 0
      if(dim(clusters[[i]])[1] > 1){
        a = mean(disim(matrix(clusters[[i]][-j,],ncol=2),point))
      }
      if(length(clusters[-i])>0){
        b = mean(disim(matrix(clusters[-i][[1]],ncol=2),point))
        if(length(clusters[-i])>1){
          for (k in 2:(length(clusters[-i]))) {
            bb = mean(disim(matrix(clusters[-i][[k]],ncol=2),point))
            if(bb<b){
              b = bb
            }
          }
        }
      }
      
      if(a<b){
        s = s + 1 - a/b
      }
      else{
        s = s + b/a - 1
      }
      iter = iter +1
    }
  }
  s = s/iter
  return(s)
}


#-------------------------------------------------------------------------------#
# 		Fontion de calcul de l'inertie intra classe d'une partition               #
#-------------------------------------------------------------------------------#
# Entrée : - clusters, partition dont on souhaite calculé le inertie       	    #
#	   - disimilarity, fonction de generation de la matrice de disimilarite	      #
#	   - up_centroid, fonction de calcul des centroides			                      #
# Sortie : - l'inertie intra classe de la partition	                            #
#-------------------------------------------------------------------------------#
inertie_intra <- function(clusters,disimilarity,up_centroid){
  inertie_intra = 0
  for (i in 1:length(clusters)) { 
    nb_var = dim(clusters[[i]])[2]
    cluster = matrix(clusters[[i]],ncol=nb_var)
    centroid = matrix(up_centroid(cluster),ncol=2)
    disim = disimilarity(cluster,centroid)
    inertie = sum(disim)/length(disim)
    inertie_intra = inertie_intra + inertie
  }
  return(inertie_intra)
}


#-------------------------------------------------------------------------------#
#  Fontion de generation de la matrice de disimilarite pour les données du CSV  #
#-------------------------------------------------------------------------------#
# Entrées : - d, la matrice des données						                              #
#	    - centroids, la matrices des centroides                           				#
# Sortie  : - mat, la matrice des disimilarite avec sur l'axe des 		          #
#	            X les centroides et sur l'axe des Y les objets 			              #
#-------------------------------------------------------------------------------#
disimilarite <- function(d,centroids){
  nb_centroid = dim(centroids)[1]
  ones = t(rep(1,nb_centroid))
  x = c(t(d[,1] %*% ones))
  y = c(t(d[,2] %*% ones))
  d_euclid = (x - centroids[,1])^2 + (y - centroids[,2])^2
  mat=matrix(d_euclid,ncol=nb_centroid,byrow=TRUE)
  return(mat)
}


#---------------------------------------------------------------#
# Fontion de mise à jour d'un centroide pour les données du CSV	#
#---------------------------------------------------------------#
# Entrée : - cluster, une classe				                        #
# Sortie : - centroid le nouveau centroid de cette classe 	    #
#---------------------------------------------------------------#
update_centroid <- function(cluster){
  x = mean(cluster[,1])
  y = mean(cluster[,2])
  centroid = c(x,y)
  return(centroid)
}


#-------------------------------------------------------------------------------#
# 		Fontion de calcul du SSE de 2 classes 				                            #
#-------------------------------------------------------------------------------#
# Entrée : - clusters, liste des 2 classes dont on souhaite calculé le SSE 	    #
#	   - disimilarity, fonction degeneration de la matrice de disimilarite	      #
#	   - up_centroid, fonction de calcul des centroides			                      #
# Sortie : - la somme des SSE des 2 classes				 	                            #
#-------------------------------------------------------------------------------#
sse_two_clusters <- function(clusters,disimilarity,up_centroid){
  nb_var = dim(data)[2]
  cluster1 = matrix(clusters[[1]],ncol=nb_var)
  centroid1 = matrix(up_centroid(cluster1),ncol=2)
  disim1 = disimilarity(cluster1,centroid1)
  
  if(length(clusters)>1){
    cluster2 =  matrix(clusters[[2]],ncol=nb_var)
    centroid2 = matrix(up_centroid(cluster2),ncol=2)
    disim2 = disimilarity(cluster2,centroid2)
  }
  else{
    disim2 = 0
  }
  
  return(sum(disim1)+sum(disim2))
}


#-------------------------------------------------------------------------------#
# 		Fontion d'affichage d'une partition    				                            #
#-------------------------------------------------------------------------------#
# Entrée : - clusters, la partition a affiché                              	    #
#	   - title, le titre du graphe                        	                      #
#-------------------------------------------------------------------------------#
plot_clusters <- function(clusters,title){
  plot(0,0,main=title, xlab="X", ylab="Y",xlim=c(0,10), ylim=c(0,10),col=0)
  for (i in 1:length(clusters)) {
    points(clusters[[i]][,1],clusters[[i]][,2],pch=i,col=i+1)
  }
}



"******************************
*     Programme principal     *
*******************************"

# Chargement des données depuis un fichier CSV
data = data.matrix(read.csv("data.csv")) 


writeLines("  ================")
writeLines("  Tests du K-Means")
writeLines("  ================")
result = kmeans_direct(data,disimilarite,update_centroid,3)
print(result)
plot_clusters(result$clusters,"Partition K-Means")


writeLines("  ===================================================")
writeLines("  Tests du K-Means avec choix des centroides initiaux")
writeLines("  ===================================================")
c = matrix(c(data[1,],data[4,],data[7,]),ncol=2,byrow=TRUE) # Les centroides sont les points 1, 4 et 7 du jeu de données 
result = kmeans_direct(data,disimilarite,update_centroid,3,centroids=c)
print(result)
plot_clusters(result$clusters,"Partition K-Means avec choix des centroides initiaux")


writeLines("  ============================")
writeLines("  Tests du K-Means Incremental")
writeLines("  ============================")
result = kmeans_incremental(data,disimilarite,update_centroid,3)
print(result)
plot_clusters(result$clusters,"Partition K-Means incremental")


writeLines("  ===========================")
writeLines("   Tests du Bisecting K-Means")
writeLines("  ===========================")
result = bisecting_kmeans(data,disimilarite,update_centroid,sse_two_clusters,3,nb_iter=50)
print(result)
plot_clusters(result$clusters,"Partition Bisecting K-Means")


writeLines("  ===========================")
writeLines("  Tests du MST based divisive")
writeLines("  ===========================")
result = mst_based_divisive(data,disimilarite,update_centroid,3)
print(result)
plot_clusters(result$clusters,"Partition MST based divisive")


writeLines("  ==========================")
writeLines("  Tests du Fuzzy K-Means m=2")
writeLines("  ==========================")
result = fuzzy_kmeans(data,2,disimilarite,0.001,3,maxit=20)
print(result)

writeLines("  ==========================")
writeLines("  Tests du Fuzzy K-Means m=3")
writeLines("  ==========================")
result = fuzzy_kmeans(data,4,disimilarite,0.001,3,maxit=20)
print(result)
