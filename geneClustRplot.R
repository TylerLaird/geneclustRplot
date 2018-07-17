##setting list.of.colors if desired----
library(randomcoloR)
set_colors<-function(file,seed.value=1){
  file<-read.table(file,
                   sep='\t',header=T,fill=T,quote='"')
  set.seed(seed.value)
  colors<-as.data.frame(t(distinctColorPalette(k=length(unique(file$Gene.Product.Name)))))
  colnames(colors)<-unique(file$Gene.Product.Name)
  colors[1,]<-apply(colors,2, function(x) as.character(x))
  return(colors)}

##gene arrow----
Gene_Arrow<-function(input,filetype='IMG',close.gene.setback=150,seed.value=1, list.of.colors=NULL){
  if(filetype=='FASTA'){
    fasta_file<-paste(readLines(input),collapse = '\n')
    headers<-unlist(str_extract_all(fasta_file,'>.+'))
    Scaffold.Name<-input
    Gene.Product.Name<-str_extract(headers,'(?<=protein=).[^]]+(?=])')
    Start.Coord<-as.numeric(str_extract(headers,'(?<=location=|location=complement\\()\\d+(?=\\.)'))
    End.Coord<-as.numeric(str_extract(headers,'(?<=\\.\\.|\\>)\\d+(?=]|)') )
    Scaffold.Length..bp.<-max(End.Coord)
    Strand<- ifelse(grepl('location=complement\\(',headers), '-','+')
    input<-cbind.data.frame(Scaffold.Name,Gene.Product.Name,Start.Coord,End.Coord,Scaffold.Length..bp.,Strand)
  }  
   input<-input[order(input$Start.Coord),]
  plot(NA, xlim=c(1,max(as.numeric(input$End.Coord))), ylim=c(-10,200),axes=F, xlab=NA, ylab=NA)
  lines(c(1,max(as.numeric(input$Scaffold.Length..bp.))),c(0,0) )
  #lines(c(1,max(as.numeric(input$Scaffold.Length..bp.))),c(-5,-5) )
  
  max_min_diff=max(input$End.Coord)-min(input$Start.Coord)
  text(c(seq(0,max(input$End.Coord),10^ceiling(log10(max_min_diff))/10   ) ), c(-8),
       gsub('0000','0k',as.character(c(seq(0,max(input$End.Coord),10^ceiling(log10(max_min_diff))/10 )))), cex=0.8 )
  
  #text(c(seq(0,100000,10000)), c(-8), gsub('0000','0k',as.character(c(seq(0,100000,10000)))), cex=0.8 )
  #text(c(seq(0,100000,10000)), c(-6), as.character('|'), cex=0.7 )
  
  
  segments(c(seq(0,max(input$End.Coord),10^ceiling(log10(max_min_diff))/10 )),c(-5),c(seq(0,max(input$End.Coord),10^ceiling(log10(max_min_diff))/10)),c(-3))
  title(input$Scaffold.Name[1],cex.main=1)
  ###
  set.seed(seed.value)
  if(is.null(list.of.colors)){
    colors<-as.data.frame(t(distinctColorPalette(k=length(unique(input$Gene.Product.Name)))))
    colnames(colors)<-unique(input$Gene.Product.Name)
    colors[1,]<-apply(colors,2, function(x) as.character(x)) }
  else{colors<-list.of.colors}
  #plotting labels taking into account closely spaced labels
  for(i in seq(1:nrow(input))){
    len<-input[i,]$End.Coord-input[i,]$Start.Coord
    mean_point<-((as.numeric(input[i,]$Start.Coord)+as.numeric(input[i,]$End.Coord))/2)
    mean_point_ahead<-((as.numeric(input[i+1,]$Start.Coord)+as.numeric(input[i+1,]$End.Coord))/2)
    mean_point_behind<-((as.numeric(input[i-1,]$Start.Coord)+as.numeric(input[i-1,]$End.Coord))/2)
    if ( i==nrow(input) ||(mean_point_ahead-mean_point >= 400|| input$Scaffold.Length..bp.<=30000)){
      segments(c(mean_point),c(0),c(mean_point),c(8))
      text(c(mean_point),c(9),input[i,]$Gene.Product.Name,cex=0.63,pos=4,srt=90,offset = 0)
    }
    else{
      segments(c(mean_point),c(0),c(mean_point),c(5))
      segments(c(mean_point),c(5),c(mean_point-close.gene.setback),c(8))
      text(c(mean_point-close.gene.setback),c(9),input[i,]$Gene.Product.Name,cex=0.63,pos=4,srt=90,offset = 0)
    }
  }
  ### plotting polygons
  apply(input, 1, function(x){ 
    mean_point<-((as.numeric(x['Start.Coord'])+as.numeric(x['End.Coord']))/2)
    span<-as.numeric(x['End.Coord'])-as.numeric(x['Start.Coord'])
    if(x['Strand']=='-'){
      polygon(c(as.numeric(x['Start.Coord'])+(span*0.2),x['Start.Coord'],as.numeric(x['Start.Coord'])+(span*0.2),x['End.Coord'],x['End.Coord'],as.numeric(x['Start.Coord'])+(span*0.2)),
              c(-3,0,3,3,-3,-3),col=as.character(colors[1,x['Gene.Product.Name']]))}
    else{
      polygon(c(x['Start.Coord'],x['Start.Coord'],as.numeric(x['End.Coord'])-(span*0.2),x['End.Coord'],as.numeric(x['End.Coord'])-(span*0.2),x['Start.Coord']),
              c(-3,3,3,0,-3,-3),col=as.character(colors[1,x['Gene.Product.Name']]))}
  })
}

