#!/usr/bin/Rscript

# $1 == $wd ; $2 == $inv

library(ggplot2)

inputs <- commandArgs(T)

inputs <- c("C:/Users/Rikip/OneDrive/Escritorio/TestsStd", "HsInv0370")

setwd(paste(inputs, collapse = "/"))

## Read all files

tc <- read.table(file = "Resultados/TablaCandidatos", sep = "\t", header = F)
colnames(tc) <- c("Read", "A", "B", "C", "D")
rownames(tc) <- tc$Read

mapinfo <- read.table(file = "Resultados/MapInfo", sep = "\t", header = F)
colnames(mapinfo) <- c("Read", "Sonda", "Signo")
mapinfo$Signo <- as.character(mapinfo$Signo)

sonErr <- unique(read.table("Resultados/sondaErronea", sep = "\t", header = F)$V1)

dupReads <- rownames(tc)
'%!in%' <- Negate('%in%')
dupReads <- dupReads[dupReads %!in% sonErr]

tablaSignos <- tc

## Leer info sobre sondas y distancias

sondas <- new.env()
for(sonda in c("A", "B", "C", "D")){
  sondas[[sonda]] <- read.table(file = paste0("CoordenadasRelativas/",sonda), header = F, sep = "\t")
}

DistSondas <- new.env()
relDir <- "Resultados/Distancias/Reads"

for (read in dupReads){
  sonPres <- dir(paste(relDir, read, "BLASTCoord" ,sep = "/"))
  if(length(sonPres) > 1){
    for(sp in sonPres){
      DistSondas[[read]][[sp]] <- read.table(file = paste(relDir, read, "BLASTCoord", sp, sep = "/"), 
                                             sep = "\t", header = F, stringsAsFactors = F)
    }
  }
}

seqInfo <- read.table(file = "seqInfo", sep = "\t", header = F)

dupReads <- ls(DistSondas)

tablaSignos <- tablaSignos[dupReads,]
tc <- tc[dupReads,]

mapinfo <- mapinfo[mapinfo$Read %in% dupReads,]

## Filtro los BC

nd <- c()

for(read in dupReads){
  sMap <- ls(DistSondas[[read]])
  
  chA <- "A" %in% sMap
  chD <- "D" %in% sMap
  
  if(sum(c(chA, chD)) > 0){
    nd <- c(nd, read)
  }
}

DistSondas <- mget(nd, envir = DistSondas)

dupReads <- ls(DistSondas)

tc <- tc[dupReads,]
mapinfo <- mapinfo[mapinfo$Read %in% dupReads,]

## Itero por cada read para ver cuantos signos tiene y cómo son.

for(sonda in c("A", "B", "C", "D")){
  tablaSignos[,paste0("signo", sonda)] <- NA
  
  for(read in dupReads){
    
    if(tc[read,sonda]){
      sInfo <- mapinfo[mapinfo$Read == read & mapinfo$Sonda == sonda, ]
      
      if(dim(sInfo)[1] == 1){
        tablaSignos[read, paste0("signo", sonda)] <- sInfo$Signo[1]
      }
      else{
        nsig <- dim(sInfo)[1]
        spos <- sum(sInfo$Signo == "+")
        
        if(spos/nsig == 1){
          tablaSignos[read, paste0("signo", sonda)] <- "+"
        }
        else if(spos/nsig == 0){
          tablaSignos[read, paste0("signo", sonda)] <- "-" 
        }
        else{        
          if(spos/nsig < 0.5){
            tablaSignos[read, paste0("signo", sonda)] <- "P-"
          } else if (spos/nsig > 0.5){
            tablaSignos[read, paste0("signo", sonda)] <- "P+"
          } else {
            tablaSignos[read, paste0("signo", sonda)] <- "NC"
          }
        }
      }
      
    }
  }
}

tablaSignos[is.na(tablaSignos)] <- "0"

## Representar Plot

refPlot <- as.data.frame(t(unlist(c(1, as.numeric(seqInfo[3]-seqInfo[2]), 
             sondas[["A"]][1], sondas[["A"]][2],
             sondas[["B"]][1], sondas[["B"]][2],
             sondas[["C"]][1], sondas[["C"]][2],
             sondas[["D"]][1], sondas[["D"]][2]))))

colnames(refPlot) <- c("Init", "End", "Ai", "Ae", "Bi", "Be", "Ci", "Ce", "Di", "De")

dataPlot <- as.data.frame(dupReads, stringsAsFactors = F)
rownames(dataPlot) <- dataPlot$dupReads
dataPlot$ReadLength <- NA

for(sonda in c("A", "B", "C", "D")){
  dataPlot[,paste0(sonda, "i")] <- NA
  dataPlot[,paste0(sonda, "e")] <- NA
}

dataPlot$AB <- NA
dataPlot$AC <- NA
dataPlot$AD <- NA
dataPlot$BD <- NA
dataPlot$CD <- NA
dataPlot$ABr <- NA
dataPlot$ACr <- NA
dataPlot$ADr <- NA
dataPlot$BDr <- NA
dataPlot$CDr <- NA

for(read in dupReads){
  sonPres <- ls(DistSondas[[read]])
  dataPlot[read, "ReadLength"] <- DistSondas[[read]][[sonPres[1]]]$V3
  
  for(sonda in c("A", "B", "C", "D")){
    if(sonda %in% sonPres){
      dataPlot[read, paste0(sonda, "i")] <- DistSondas[[read]][[sonda]]$V1
      dataPlot[read, paste0(sonda, "e")] <- DistSondas[[read]][[sonda]]$V2
    }
  }
  
  if(!is.na(dataPlot[read,"Ai"])){
    if(!is.na(dataPlot[read,"Bi"])){
      dataPlot[read, "AB"] <- ((dataPlot[read, "Be"] + dataPlot[read, "Bi"])/2) -
        ((dataPlot[read, "Ae"] + dataPlot[read, "Ai"])/2)
      dataPlot[read, "ABr"] <- abs(dataPlot[read, "AB"])/
        (((refPlot$Be + refPlot$Bi)/2)-((refPlot$Ae + refPlot$Ai)/2))
    }
    if(!is.na(dataPlot[read,"Ci"])){
      dataPlot[read, "AC"] <- ((dataPlot[read, "Ce"] + dataPlot[read, "Ci"])/2) -
        ((dataPlot[read, "Ae"] + dataPlot[read, "Ai"])/2)
      dataPlot[read, "ACr"] <- abs(dataPlot[read, "AC"])/
        (((refPlot$Ce + refPlot$Ci)/2)-((refPlot$Ae + refPlot$Ai)/2))
    }
    if(!is.na(dataPlot[read,"Di"])){
      dataPlot[read, "AD"] <- ((dataPlot[read, "De"] + dataPlot[read, "Di"])/2) -
        ((dataPlot[read, "Ae"] + dataPlot[read, "Ai"])/2)
      dataPlot[read, "ADr"] <- abs(dataPlot[read, "AD"])/
        (((refPlot$De + refPlot$Di)/2)-((refPlot$Ae + refPlot$Ai)/2))
    }
  }
  if(!is.na(dataPlot[read,"Di"])){
    if(!is.na(dataPlot[read,"Bi"])){
      dataPlot[read, "BD"] <- ((dataPlot[read, "De"] + dataPlot[read, "Di"])/2) -
        ((dataPlot[read, "Be"] + dataPlot[read, "Bi"])/2)
      dataPlot[read, "BDr"] <- abs(dataPlot[read, "BD"])/
        (((refPlot$De + refPlot$Di)/2)-((refPlot$Be + refPlot$Bi)/2))
    }
    if(!is.na(dataPlot[read,"Ci"])){
      dataPlot[read, "CD"] <- ((dataPlot[read, "De"] + dataPlot[read, "Di"])/2) -
        ((dataPlot[read, "Ce"] + dataPlot[read, "Ci"])/2)
      dataPlot[read, "CDr"] <- abs(dataPlot[read, "CD"])/
        (((refPlot$De + refPlot$Di)/2)-((refPlot$Ce + refPlot$Ci)/2))
    }
  }
}

repGraf <- data.frame(matrix(nrow = 0, ncol = 13))

for(read in dupReads){
  forward <- sum(dataPlot[read, c("AB", "AC", "AD", "BD", "CD")], na.rm = T) > 0
  rRl <- dataPlot[read, "ReadLength"]
  
  if(!is.na(dataPlot[read, "Ai"])){
    rAi <- refPlot$Ai
    rAe <- refPlot$Ae
    
    if(!is.na(dataPlot[read, "Bi"])){
      rBi <- rAi + abs(dataPlot[read, "AB"])
      rBe <- rAe + abs(dataPlot[read, "AB"])
    } else {
      rBi <- NA
      rBe <- NA
    }
    if(!is.na(dataPlot[read, "Ci"])){
      rCi <- rAi + abs(dataPlot[read, "AC"])
      rCe <- rAe + abs(dataPlot[read, "AC"])
    } else {
      rCi <- NA
      rCe <- NA
    }
    if(!is.na(dataPlot[read, "Di"])){
      rDi <- rAi + abs(dataPlot[read, "AD"])
      rDe <- rAe + abs(dataPlot[read, "AD"])
    } else {
      rDi <- NA
      rDe <- NA
    }
    
    if(forward){
      rRi <- refPlot$Ai - dataPlot[read, "Ai"]
      if(rRi < 0){
        rRi <- 0
      }
      rRe <- refPlot$Ae + (rRl - dataPlot[read, "Ae"])
      if(rRe > refPlot$End){
        rRe <- refPlot$End
      }
    } else {
      rRi <- refPlot$Ai - (rRl - dataPlot[read, "Ai"])
      if(rRi < 0){
        rRi <- 0
      }
      rRe <- refPlot$Ae + dataPlot[read, "Ae"]
      if(rRe > refPlot$End){
        rRe <-  refPlot$End
      }
    }
    
  } else {
    rAi <- NA
    rAe <- NA
    
    rDi <- refPlot$Di
    rDe <- refPlot$De
    
    if(!is.na(dataPlot[read, "Bi"])){
      rBi <- rDi - abs(dataPlot[read, "BD"])
      rBe <- rDe - abs(dataPlot[read, "BD"])
    } else {
      rBi <- NA
      rBe <- NA
    }
    if(!is.na(dataPlot[read, "Ci"])){
      rCi <- rDi - abs(dataPlot[read, "CD"])
      rCe <- rDe - abs(dataPlot[read, "CD"])
    } else {
      rCi <- NA
      rCe <- NA
    }
    
    if(forward){
      rRi <- refPlot$Di - dataPlot[read, "Di"]
      if(rRi < 0){
        rRi <-  0
      }
      rRe <- refPlot$De + (rRl - dataPlot[read, "De"])
      if(rRe > refPlot$End){
        rRe <-  refPlot$End
      }
    } else {
      rRi <- refPlot$Di - (rRl - dataPlot[read, "Di"])
      if(rRi < 0){
        rRi <-  0
      }
      rRe <- refPlot$De + dataPlot[read, "De"]
      if(rRe > refPlot$End){
        rRe <-  refPlot$End
      }
    }
  }
  
  Or <- NA
  repGraf <- rbind(repGraf, c(read, rAi, rAe, rBi, rBe, rCi, rCe, rDi, rDe, rRi, rRe, rRl, Or))
  
}

repGraf <- rbind(repGraf, c("Reference", refPlot$Ai, refPlot$Ae, refPlot$Bi, refPlot$Be, refPlot$Ci, refPlot$Ce,
                            refPlot$Di, refPlot$De, NA, NA, NA, NA))
colnames(repGraf) <- c("Read", "Ai", "Ae", "Bi", "Be", "Ci", "Ce", "Di", "De", "Ri", "Re", "Rl", "Or")
rownames(repGraf) <- repGraf$Read

for(cn in colnames(repGraf)){
  repGraf[,cn] <- as.numeric(repGraf[,cn])
}

rG <- repGraf
rG[is.na(rG)] <- -40000

hGraph <- 200 + 100*dim(repGraf)[1]
iterGraf <- c("Reference", dupReads)

png(filename = "ClevelandDistancias.png", height = hGraph, width = 2000)
  ggplot(repGraf) +
    # geom_point(aes(x = (Ae+Ai)/2, y = 1), color = "#90af43", size = 10) +
    # geom_point(aes(x = (Be+Bi)/2, y = 1), color = "#932e1f", size = 10) +
    # geom_point(aes(x = (Ce+Ci)/2, y = 1), color = "#e4a21a", size = 10) +
    # geom_point(aes(x = (De+Di)/2, y = 1), color = "#424daa", size = 10) +
    geom_rect(aes(xmin = 0, xmax = refPlot[1,2], ymin = dim(repGraf)[1]-0.05, ymax = dim(repGraf)[1]+0.05), fill = "#222222", color = NA, size = 1) +
    geom_rect(aes(xmin = Ri, xmax = Re, ymin = seq(length(iterGraf))-0.05, ymax = seq(length(iterGraf))+0.05), fill = "#d3d3d3", color = NA, size = 1) +
    geom_rect(aes(xmin = Ai, xmax = Ae, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#90af43", color = NA, size = 4) +
    geom_rect(aes(xmin = Bi, xmax = Be, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#932e1f", color = NA, size = 4) +
    geom_rect(aes(xmin = Ci, xmax = Ce, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#e4a21a", color = NA, size = 4) +
    geom_rect(aes(xmin = Di, xmax = De, ymin = seq(length(iterGraf))-0.3, ymax = seq(length(iterGraf))+0.3), fill = "#424daa", color = NA, size = 4) +
    theme_classic() + 
    xlim(0,refPlot[1,2]) + 
    ylim(0,dim(repGraf)[1]+0.5)

dev.off()





