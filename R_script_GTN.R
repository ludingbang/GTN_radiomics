#Heatmap (Figure 2)
plot(hclust(as.dist (1 - cor (t(scale(GTN_radiomics_clinical_all_doppler[,2:165])) / 2)),method="complete"),label=GTN_radiomics_clinical_all_doppler$Class,main=c("Dendrogram of all cases","N=221, hcluster by Pearson correlation"))
plot(hclust(dist (((scale(GTN_radiomics_clinical_all_doppler[,2:165])))),method="complete"),label=GTN_radiomics_clinical_all_doppler$Class,main=c("Dendrogram of serous cases","N=248, hcluster by Euclidean distance"))

my_palette <- colorRampPalette(c("red", "yellow","blue"))(n = 504)
cols1 = brewer.pal(4, "Blues")
pal1 = colorRampPalette(cols1)
LeftPIorder = findInterval(GTN_radiomics_clinical_all_doppler$Left.PI,sort(GTN_radiomics_clinical_all_doppler$Left.PI))
cols2 = brewer.pal(4, "Greens")
pal2 = colorRampPalette(cols2)
RightPIorder = findInterval(GTN_radiomics_clinical_all_doppler$Right.PI,sort(GTN_radiomics_clinical_all_doppler$Right.PI))

cols3 = brewer.pal(4, "Reds")
pal3 = colorRampPalette(cols3)
WHOhCGScoreneworder = findInterval(GTN_radiomics_clinical_all_doppler$WHOhCGScorenew,sort(GTN_radiomics_clinical_all_doppler$WHOhCGScorenew))

heatmap.plus(t(scale(GTN_radiomics_clinical_all_doppler[,2:165],scale=T)), cexCol=0.1,cexRow=0.5,distfun=function(c) as.dist (1 - cor (t(c) / 2)),col=my_palette,labRow=FALSE,labCol=T,breaks=c(seq(-15,-0.8,length=101),seq(-0.8,-0.1,length=101),seq(-0.1,0.1,length=101),seq(0.1,0.8,length=101),seq(0.8,15,length=101)),ColSideColors=cbind(GTN_radiomics_clinical_all_doppler$Clustercolor,GTN_radiomics_clinical_all_doppler$Classcolor,GTN_radiomics_clinical_all_doppler$FIGOcolorbi,pal1(ncol(GTN_radiomics_clinical_all_doppler))[LeftPIorder],pal2(ncol(GTN_radiomics_clinical_all_doppler))[RightPIorder],pal3(ncol(GTN_radiomics_clinical_all_doppler))[WHOhCGScoreneworder]),RowSideColors=cbind(dopplerradiomics$V2,dopplerradiomics$V2))

#Volcano plot
set.seed(9) 
train1=sample(seq_len(nrow(x_doppler_scale)),size = 86)

x_training=x_doppler_scale[train1,]

test=sapply(165:328,function(n)t.test(x_training[,ID2[n]]~x_training[,"class_doppler"])$p.value)

test=as.data.frame(as.matrix(test))
row.names(test)=names(x_training)[1:164]
test$fdr=p.adjust(test[,1],method="fdr")
names(test)[1]="p"

test$fc=sapply(165:328,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$estimate[2])-sapply(165:328,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$estimate[1])

n=cbind(intercept=1,group=x_training[,"class_doppler"])

test=topTable(eBayes(lmFit(t(x_training[,ID2[165:328]]),design= n)),coef=2,n=164)

ggplot(test, aes(logFC, -log(P.Value,10),color=logFC)) +
  geom_point(shape = 16, size = 4, show.legend = FALSE, alpha = .6) +
  theme_minimal() + geom_hline(yintercept=-log(0.01,10),linetype=2,color="grey55")+geom_vline(xintercept=0,linetype=2,color="grey55")+geom_text( mapping=aes(label=ID), size=2, vjust=-1, hjust=0.5)+xlim(-1.5,1.5)+ylim(0, 5)


#Both ultrasound+Doppler 100 times random split of dataset (Supplementary Figure 1)
GTN_ren_sumamry_twosets=matrix(0, ncol = 8, nrow = 100)
GTN_ren_sumamry_twosets=as.data.frame(GTN_ren_sumamry_twosets)
names(GTN_ren_sumamry_twosets)=c("p005_num","fdr25_num","auc_train","auc_val1","auc_val2")

for (i in 1:100){
  set.seed(i) 
  train1=sample(seq_len(nrow(x_ultrasound_doppler_scale)),size = 74)
  
  x_training=x_ultrasound_doppler_scale[train1,]
  
  val1=sample(seq_len(nrow(x_ultrasound_doppler_scale[-train1,])),size = 74)
  x_val1=x_ultrasound_doppler_scale[-train1,][val1,]
  x_val2=x_ultrasound_doppler_scale[-train1,][-val1,]
  
  ttest_ren_training=sapply(1:328,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$p.value)
  
  ttest_ren_training=as.data.frame(as.matrix(ttest_ren_training))
  row.names(ttest_ren_training)=names(x_training)[2:329]
  ttest_ren_training$fdr=p.adjust(ttest_ren_training[,1],method="fdr")
  names(ttest_ren_training)[1]="p"
  
  if (dim(subset(ttest_ren_training,p<0.05))[1]<2){
    GTN_ren_sumamry_twosets[i,1]="No significant univariable features"
    GTN_ren_sumamry_twosets[i,2]="No significant univariable features"
  } else {
    GTN_ren_sumamry_twosets[i,1]=dim(subset(ttest_ren_training,p<0.05))[1]
    GTN_ren_sumamry_twosets[i,2]=dim(subset(ttest_ren_training,fdr<0.25))[1]
    
    cvfit_NRvsR = cv.glmnet(x=data.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), y=as.factor(x_training$class), family='binomial', alpha=1)
    
    GTN_lasso_val1_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val1[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_val2_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val2[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_training_ren=predict(cvfit_NRvsR, newx = as.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    
    if (median(GTN_lasso_val1_ren)==0){
      GTN_ren_sumamry_twosets[i,3]="No lambda.min"
      GTN_ren_sumamry_twosets[i,4]="No lambda.min"
      GTN_ren_sumamry_twosets[i,5]="No lambda.min"
    } else {
      
      pr <- prediction(GTN_lasso_val1_ren, x_val1$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_twosets[i,4]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_val2_ren, x_val2$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_twosets[i,5]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_training_ren, x_training$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_twosets[i,3]=auc@y.values[[1]]
    }
  }
}


#Ultrasound (overlapping cases with doppler) 100times random split of dataset (Supplementary Figure 1)
GTN_ren_sumamry_set1=matrix(0, ncol = 8, nrow = 100)
GTN_ren_sumamry_set1=as.data.frame(GTN_ren_sumamry_set1)
names(GTN_ren_sumamry_set1)=c("p005_num","fdr25_num","auc_train","auc_val1","auc_val2")

for (i in 1:100){
  set.seed(i) 
  train1=sample(seq_len(nrow(x_ultrasound_doppler_scale)),size = 74)
  
  x_training=x_ultrasound_doppler_scale[train1,]
  
  val1=sample(seq_len(nrow(x_ultrasound_doppler_scale[-train1,])),size = 74)
  x_val1=x_ultrasound_doppler_scale[-train1,][val1,]
  x_val2=x_ultrasound_doppler_scale[-train1,][-val1,]
  
  ttest_ren_training=sapply(1:164,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$p.value)
  
  ttest_ren_training=as.data.frame(as.matrix(ttest_ren_training))
  row.names(ttest_ren_training)=names(x_training)[2:165]
  ttest_ren_training$fdr=p.adjust(ttest_ren_training[,1],method="fdr")
  names(ttest_ren_training)[1]="p"
  
  if (dim(subset(ttest_ren_training,p<0.05))[1]<2){
    GTN_ren_sumamry_set1[i,1]="No significant univariable features"
    GTN_ren_sumamry_set1[i,2]="No significant univariable features"
  } else {
    GTN_ren_sumamry_set1[i,1]=dim(subset(ttest_ren_training,p<0.05))[1]
    GTN_ren_sumamry_set1[i,2]=dim(subset(ttest_ren_training,fdr<0.25))[1]
    
    cvfit_NRvsR = cv.glmnet(x=data.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), y=as.factor(x_training$class), family='binomial', alpha=1)
    
    GTN_lasso_val1_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val1[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_val2_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val2[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_training_ren=predict(cvfit_NRvsR, newx = as.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    
    if (median(GTN_lasso_val1_ren)==0){
      GTN_ren_sumamry_set1[i,3]="No lambda.min"
      GTN_ren_sumamry_set1[i,4]="No lambda.min"
      GTN_ren_sumamry_set1[i,5]="No lambda.min"
    } else {
      
      pr <- prediction(GTN_lasso_val1_ren, x_val1$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set1[i,4]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_val2_ren, x_val2$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set1[i,5]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_training_ren, x_training$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set1[i,3]=auc@y.values[[1]]
    }
  }
}

#doppler (overlapping cases with ultrasound) 100times random split of dataset (Supplementary Figure 1)
GTN_ren_sumamry_set2=matrix(0, ncol = 8, nrow = 100)
GTN_ren_sumamry_set2=as.data.frame(GTN_ren_sumamry_set2)
names(GTN_ren_sumamry_set2)=c("p005_num","fdr25_num","auc_train","auc_val1","auc_val2")

for (i in 1:100){
  set.seed(i) 
  train1=sample(seq_len(nrow(x_ultrasound_doppler_scale)),size = 74)
  
  x_training=x_ultrasound_doppler_scale[train1,]
  
  val1=sample(seq_len(nrow(x_ultrasound_doppler_scale[-train1,])),size = 74)
  x_val1=x_ultrasound_doppler_scale[-train1,][val1,]
  x_val2=x_ultrasound_doppler_scale[-train1,][-val1,]
  
  ttest_ren_training=sapply(165:328,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$p.value)
  
  ttest_ren_training=as.data.frame(as.matrix(ttest_ren_training))
  row.names(ttest_ren_training)=names(x_training)[166:329]
  ttest_ren_training$fdr=p.adjust(ttest_ren_training[,1],method="fdr")
  names(ttest_ren_training)[1]="p"
  
  if (dim(subset(ttest_ren_training,p<0.05))[1]<2){
    GTN_ren_sumamry_set2[i,1]="No significant univariable features"
    GTN_ren_sumamry_set2[i,2]="No significant univariable features"
  } else {
    GTN_ren_sumamry_set2[i,1]=dim(subset(ttest_ren_training,p<0.05))[1]
    GTN_ren_sumamry_set2[i,2]=dim(subset(ttest_ren_training,fdr<0.25))[1]
    
    cvfit_NRvsR = cv.glmnet(x=data.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), y=as.factor(x_training$class), family='binomial', alpha=1)
    
    GTN_lasso_val1_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val1[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_val2_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val2[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_training_ren=predict(cvfit_NRvsR, newx = as.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    
    if (median(GTN_lasso_val1_ren)==0){
      GTN_ren_sumamry_set2[i,3]="No lambda.min"
      GTN_ren_sumamry_set2[i,4]="No lambda.min"
      GTN_ren_sumamry_set2[i,5]="No lambda.min"
    } else {
      
      pr <- prediction(GTN_lasso_val1_ren, x_val1$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set2[i,4]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_val2_ren, x_val2$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set2[i,5]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_training_ren, x_training$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set2[i,3]=auc@y.values[[1]]
    }
  }
}


#doppler all cases 100 times random split
GTN_ren_sumamry_set2_all=matrix(0, ncol = 8, nrow = 100)
GTN_ren_sumamry_set2_all=as.data.frame(GTN_ren_sumamry_set2_all)
names(GTN_ren_sumamry_set2_all)=c("p005_num","fdr25_num","auc_train","auc_val1","auc_val2")

for (i in 1:100){
  set.seed(i) 
  train1=sample(seq_len(nrow(x_doppler_scale)),size = 86)
  
  x_training=x_doppler_scale[train1,]
  
  val1=sample(seq_len(nrow(x_doppler_scale[-train1,])),size = 86)
  x_val1=x_doppler_scale[-train1,][val1,]
  x_val2=x_doppler_scale[-train1,][-val1,]
  
  ttest_ren_training=sapply(165:328,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$p.value)
  
  ttest_ren_training=as.data.frame(as.matrix(ttest_ren_training))
  row.names(ttest_ren_training)=names(x_training)[1:164]
  ttest_ren_training$fdr=p.adjust(ttest_ren_training[,1],method="fdr")
  names(ttest_ren_training)[1]="p"
  
  if (dim(subset(ttest_ren_training,p<0.05))[1]<2){
    GTN_ren_sumamry_set2_all[i,1]="No significant univariable features"
    GTN_ren_sumamry_set2_all[i,2]="No significant univariable features"
  } else {
    GTN_ren_sumamry_set2_all[i,1]=dim(subset(ttest_ren_training,p<0.05))[1]
    GTN_ren_sumamry_set2_all[i,2]=dim(subset(ttest_ren_training,fdr<0.25))[1]
    
    cvfit_NRvsR = cv.glmnet(x=data.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), y=as.factor(x_training$class), family='binomial', alpha=1)
    
    GTN_lasso_val1_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val1[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_val2_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val2[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    GTN_lasso_training_ren=predict(cvfit_NRvsR, newx = as.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
    
    if (median(GTN_lasso_val1_ren)==0){
      GTN_ren_sumamry_set2_all[i,3]="No lambda.min"
      GTN_ren_sumamry_set2_all[i,4]="No lambda.min"
      GTN_ren_sumamry_set2_all[i,5]="No lambda.min"
    } else {
      
      pr <- prediction(GTN_lasso_val1_ren, x_val1$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set2_all[i,4]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_val2_ren, x_val2$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set2_all[i,5]=auc@y.values[[1]]
      
      pr <- prediction(GTN_lasso_training_ren, x_training$class)
      auc=performance(pr, measure = "auc")
      GTN_ren_sumamry_set2_all[i,3]=auc@y.values[[1]]
    }
  }
}

#Final model selection using Doppler radiomics only (Figure 3)
set.seed(9) 
train1=sample(seq_len(nrow(x_doppler_scale)),size = 86)

x_training=x_doppler_scale[train1,]

val1=sample(seq_len(nrow(x_doppler_scale[-train1,])),size = 86)
x_val1=x_doppler_scale[-train1,][val1,]
x_val2=x_doppler_scale[-train1,][-val1,]

ttest_ren_training=sapply(165:328,function(n)t.test(x_training[,ID2[n]]~ x_training[,"class_doppler"])$p.value)

ttest_ren_training=as.data.frame(as.matrix(ttest_ren_training))
row.names(ttest_ren_training)=names(x_training)[1:164]
ttest_ren_training$fdr=p.adjust(ttest_ren_training[,1],method="fdr")
names(ttest_ren_training)[1]="p"

cvfit_NRvsR = cv.glmnet(x=data.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), y=as.factor(x_training$class), family='binomial', alpha=1)

GTN_lasso_val1_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val1[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
GTN_lasso_val2_ren=predict(cvfit_NRvsR, newx = as.matrix(x_val2[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
GTN_lasso_training_ren=predict(cvfit_NRvsR, newx = as.matrix(x_training[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")

pr_val1 <- prediction(GTN_lasso_val1_ren, x_val1$class)
prf <- performance(pr_val1, measure = "tpr", x.measure = "fpr")
plot(prf)
performance(pr_val1, measure = "auc")@y.values[[1]]

pr_val2 <- prediction(GTN_lasso_val2_ren, x_val2$class)
prf <- performance(pr_val2, measure = "tpr", x.measure = "fpr")
plot(prf)
performance(pr_val2, measure = "auc")@y.values[[1]]

pr_training <- prediction(GTN_lasso_training_ren, x_training$class)
prf <- performance(pr_training, measure = "tpr", x.measure = "fpr")
plot(prf)
performance(pr_training, measure = "auc")@y.values[[1]]

GTN_lasso_all_doppler=rbind(GTN_lasso_training_ren, GTN_lasso_val1_ren, GTN_lasso_val2_ren)
GTN_lasso_all_doppler=as.data.frame(GTN_lasso_all_doppler)

names(GTN_lasso_all_doppler)="LASSO_doppler"
GTN_lasso_all_doppler$Set=c(rep("Training",86),rep("Val1",86),rep("Val2",86))

GTN_radiomics_clinical_all_doppler=merge(GTN_radiomics_clinical_all, GTN_lasso_all_doppler,by.x="Row.names",by.y="row.names")

summary(glm(class~LASSO_doppler+FIGO3+ Left.PI+Right.PI,data= GTN_radiomics_clinical_all_doppler))

#AUC for 5 feature models
GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

GTN_glm_multi=glm(Class~ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]


GTN_glm_multi=glm(Class~ LASSO_doppler,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]


GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

GTN_glm_multi=glm(Class~ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]


GTN_glm_multi=glm(Class~ LASSO_doppler,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]


GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

GTN_glm_multi=glm(Class~ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]


GTN_glm_multi=glm(Class~ LASSO_doppler,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm, subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
prf2_val1 <- performance(pr_val1_glm, measure = "tpr", x.measure = "fpr")
plot(prf2_val1)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

#Choosing cutoffs (Figure 6)
cutoffs <- data.frame(cut=prf2_val1@alpha.values[[1]], fpr=prf2_val1@x.values[[1]], 
                      tpr=prf2_val1@y.values[[1]])

cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
head(subset(cutoffs, fpr < 0.25))

cutoffs <- data.frame(cut=prf1_val1@alpha.values[[1]], fpr=prf1_val1@x.values[[1]], 
                      tpr=prf1_val1@y.values[[1]])

cutoffs <- cutoffs[order(cutoffs$tpr, decreasing=TRUE),]
head(subset(cutoffs, fpr < 0.25))



#ggplot ROC (Figure 3 and Supplementary Figure 2)
library(plotROC)

GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm1=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm2=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ FIGO3,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm3=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ Right.PI+ Left.PI,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm4=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm5=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ LASSO_doppler,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm6=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

test=as.data.frame(cbind(subset(GTN_radiomics_clinical_all_doppler,Set=="Val1"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class")],GTN_lasso_val1_ren_glm1,GTN_lasso_val1_ren_glm2,GTN_lasso_val1_ren_glm3,GTN_lasso_val1_ren_glm4,GTN_lasso_val1_ren_glm5,GTN_lasso_val1_ren_glm6))

names(test)=c("V1","Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","FIGO","PI","hCG","Radiomics")
longtest <- melt_roc(test, "V1", c("Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","FIGO","PI","hCG","Radiomics"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)

longtest <- melt_roc(test, "V1", c("Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","Radiomics"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)


GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm1=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm2=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ FIGO3,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm3=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ Right.PI+ Left.PI,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm4=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm5=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ LASSO_doppler,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm6=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

test=as.data.frame(cbind(subset(GTN_radiomics_clinical_all_doppler,Set=="Val2"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class")],GTN_lasso_val1_ren_glm1,GTN_lasso_val1_ren_glm2,GTN_lasso_val1_ren_glm3,GTN_lasso_val1_ren_glm4,GTN_lasso_val1_ren_glm5,GTN_lasso_val1_ren_glm6))

names(test)=c("V1","Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","FIGO","PI","hCG","Radiomics")
longtest <- melt_roc(test, "V1", c("Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","FIGO","PI","hCG"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)

longtest <- melt_roc(test, "V1", c("Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","Radiomics"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)


GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm1=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm2=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ FIGO3,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm3=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ Right.PI+ Left.PI,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm4=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm5=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ LASSO_doppler,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm6=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

test=as.data.frame(cbind(subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class")],GTN_lasso_val1_ren_glm1,GTN_lasso_val1_ren_glm2,GTN_lasso_val1_ren_glm3,GTN_lasso_val1_ren_glm4,GTN_lasso_val1_ren_glm5,GTN_lasso_val1_ren_glm6))

names(test)=c("V1","Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","FIGO","PI","hCG","Radiomics")
longtest <- melt_roc(test, "V1", c("Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","FIGO","PI","hCG"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)

longtest <- melt_roc(test, "V1", c("Radiomics+PI+FIGO+hCG","PI+FIGO+hCG","Radiomics"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm5, subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

#ggroc comparing ultrasound and doppler (Supplementary Figure 2)
GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(test,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm1=predict(GTN_glm_multi, newdata = subset(test,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ LASSO_ultra+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(test1,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_ultra", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm2=predict(GTN_glm_multi, newdata = subset(test1,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_ultra", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

GTN_glm_multi=glm(Class~ LASSO_ultradop+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(test2,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_ultradop", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

GTN_lasso_val1_ren_glm3=predict(GTN_glm_multi, newdata = subset(test2,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_ultradop", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

test=as.data.frame(cbind(subset(test,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class")],GTN_lasso_val1_ren_glm1,GTN_lasso_val1_ren_glm2,GTN_lasso_val1_ren_glm3))

names(test)=c("V1","Doppler+PI+FIGO+hCG","Ultrasound+PI+FIGO+hCG","Ultra+Doppler+PI+FIGO+hCG")
longtest <- melt_roc(test, "V1", c("Doppler+PI+FIGO+hCG","Ultrasound+PI+FIGO+hCG","Ultra+Doppler+PI+FIGO+hCG"))

ggplot(longtest, aes(d = D, m = M,color=name)) + geom_roc(labels = FALSE)+ style_roc(theme = theme_gray)

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm1, subset(test,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm2, subset(test1,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_ultra", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]

pr_val1_glm <- prediction(GTN_lasso_val1_ren_glm3, subset(test2,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_ultradop", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)
performance(pr_val1_glm, measure = "auc")@y.values[[1]]


#Barplot for radiomics vs stage, class, PI, hCG (Figure 4)

test=subset(GTN_radiomics_clinical_all_doppler, FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Row.names","LASSO_doppler","Class", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]
x1=melt(test,id.var="Row.names",measure.var="LASSO_doppler")
x1=merge(x1,test[,c(1,2)],by.x=1,by.y=1)

#For barchart:
p1= ggplot(x1,aes(reorder(Row.names,LASSO_doppler),value,fill=variable))+geom_bar(stat= "identity")+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.key.size = unit(0.3,"cm"),legend.key.width =unit(0.15,"cm"),panel.background =  element_rect(fill = NA))+
  geom_bar(data=x1,stat= "identity")+
  scale_fill_brewer(palette="Set2")+geom_line(group=1,colour="darkolivegreen")

#For categorical:
x2=melt(test,id.var="Row.names",measure.var=c("WHOhCGScorenew", "FIGO3","Class"))
x2=merge(x2, test[,c(1,2)],by.x=1,by.y=1)


p2= ggplot(x2, aes(x=reorder(Row.names,LASSO_doppler), y=variable))+
  geom_tile(aes(fill = value),width=0.8, height= 0.8,colour="white") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.key.size = unit(0.3,"cm"),legend.key.width =unit(0.15,"cm")) +
  scale_fill_manual(values=c("lightcyan1","lightblue1","lightblue3","lightblue4","pink","tomato1","red","chocolate","chocolate4"),na.value="grey95")


x3=test[,c("Row.names","Right.PI", "Left.PI")]
x3[,2:3]=scale(x3[,2:3])
x3=melt(x3,id.var="Row.names",measure.var=c("Right.PI", "Left.PI"))
x3=merge(x3, test[,c(1,2)],by.x=1,by.y=1)

p3= ggplot(x3, aes(x=reorder(Row.names,LASSO_doppler), y=variable))+
  geom_tile(aes(fill = value),width=0.8, height= 0.8,colour="white", ) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank(),legend.key.size = unit(0.3,"cm"),legend.key.width =unit(0.15,"cm")) +
  scale_fill_distiller(palette="RdYlBu" , na.value="grey95")

library(ggpubr)

ggarrange(p1 ,p2, p3 + font("x.text", size = 10),
          ncol = 1, nrow = 3,align = "v", heights=c(0.2,0.4,0.3))

#Nomogram (Figure 5)
d=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]
d$FIGO3=as.factor(d$FIGO3)
d$Class=as.factor(d$Class)
d1= datadist(d)
options(datadist = 'd1')

fittest=lrm(as.factor(Class)~ LASSO_doppler+Right.PI+Left.PI+WHOhCGScorenew+FIGO3,data=d)

plot(nomogram(fittest, fun = function(x)plogis(x),lp=FALSE))


d=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class", "LASSO_doppler","FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]
d$FIGO3=as.factor(d$FIGO3)
d$Class=as.factor(d$Class)
d1= datadist(d)
options(datadist = 'd1')

fittest=lrm(as.factor(Class)~ Right.PI+Left.PI+WHOhCGScorenew+FIGO3,data=d)

plot(nomogram(fittest, fun = function(x)plogis(x),lp=FALSE))



#Confusion matrix (Figure 6)
GTN_glm_multi=glm(Class~ LASSO_doppler+ FIGO3+ Right.PI+ Left.PI+WHOhCGScorenew,data=subset(GTN_radiomics_clinical_all_doppler,Set=="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],family="binomial")

class_pred=predict(GTN_glm_multi, newdata = subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")],  type = "response")

# Calculate class probabilities: pred_class

pred_class <- as.factor(ifelse(class_pred > 0.5945965, "Resistant","Nonresistant"))

cmtrx <- confusionMatrix(pred_class, subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)

#Confusion matrix
draw_confusion_matrix <- function(cmtrx) {
  
  total <- sum(cmtrx$table)
  
  res <- as.numeric(cmtrx$table)
  
  # Generate color gradients. Palettes come from RColorBrewer.
  
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  
  getColor <- function (greenOrRed = "green", amount = 0) {
    
    if (amount == 0)
      
      return("#FFFFFF")
    
    palette <- greenPalette
    
    if (greenOrRed == "red")
      
      palette <- redPalette
    
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
    
  }
  
  # set the basic layout
  
  layout(matrix(c(1,1,2)))
  
  par(mar=c(2,2,2,2))
  
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  
  title('CONFUSION MATRIX', cex.main=2)
  
  # create the matrix
  
  classes = colnames(cmtrx$table)
  
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  
  text(195, 435, classes[1], cex=1.2)
  
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  
  text(295, 435, classes[2], cex=1.2)
  
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  
  text(245, 450, 'Actual', cex=1.3, font=2)
  
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  
  text(140, 400, classes[1], cex=1.2, srt=90)
  
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cmtrx results
  
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics
  
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  
  text(10, 85, names(cmtrx$byClass[1]), cex=1.2, font=2)
  
  text(10, 70, round(as.numeric(cmtrx$byClass[1]), 3), cex=1.2)
  
  text(30, 85, names(cmtrx$byClass[2]), cex=1.2, font=2)
  
  text(30, 70, round(as.numeric(cmtrx$byClass[2]), 3), cex=1.2)
  
  text(50, 85, names(cmtrx$byClass[5]), cex=1.2, font=2)
  
  text(50, 70, round(as.numeric(cmtrx$byClass[5]), 3), cex=1.2)
  
  text(70, 85, names(cmtrx$byClass[6]), cex=1.2, font=2)
  
  text(70, 70, round(as.numeric(cmtrx$byClass[6]), 3), cex=1.2)
  
  text(90, 85, names(cmtrx$byClass[7]), cex=1.2, font=2)
  
  text(90, 70, round(as.numeric(cmtrx$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information
  
  text(30, 35, names(cmtrx$overall[1]), cex=1.5, font=2)
  
  text(30, 20, round(as.numeric(cmtrx$overall[1]), 3), cex=1.4)
  
  text(70, 35, names(cmtrx$overall[2]), cex=1.5, font=2)
  
  text(70, 20, round(as.numeric(cmtrx$overall[2]), 3), cex=1.4)
  
}



draw_confusion_matrix(cmtrx)


# Calculate class probabilities: pred_class

pred_class <- as.factor(ifelse(class_pred > 0.6496973, "Resistant","Nonresistant"))

cmtrx <- confusionMatrix(pred_class, subset(GTN_radiomics_clinical_all_doppler,Set!="Training"& FIGO3<4&Right.PI<100000&Left.PI<100000)[,c("Class","LASSO_doppler", "FIGO3", "Right.PI", "Left.PI","WHOhCGScorenew")]$Class)

draw_confusion_matrix(cmtrx)

#BD test for specificity

BDtest(matrix(c(66,21,11,44),ncol=2), pr=0.4, conf.level = 0.95)

#Erosion and dilation (Supplementary Figure 4)
Lasso_E1=predict(cvfit_NRvsR, newx = as.matrix(MorphE1_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_E2=predict(cvfit_NRvsR, newx = as.matrix(MorphE2_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_E3=predict(cvfit_NRvsR, newx = as.matrix(MorphE3_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_E4=predict(cvfit_NRvsR, newx = as.matrix(MorphE4_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_D1=predict(cvfit_NRvsR, newx = as.matrix(MorphD1_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_D2=predict(cvfit_NRvsR, newx = as.matrix(MorphD2_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_D3=predict(cvfit_NRvsR, newx = as.matrix(MorphD3_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_D4=predict(cvfit_NRvsR, newx = as.matrix(MorphD4_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")
Lasso_no=predict(cvfit_NRvsR, newx = as.matrix(noMorph_c_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")


#Inter-observer (Supplementary Figure 5)
rep2021_results_scale=scale(rep2021_results_scale[,row.names(scalefactor_doppler)],scale=scalefactor_doppler$scale,center=scalefactor_doppler$center)

Lasso_rep=predict(cvfit_NRvsR, newx = as.matrix(rep2021_results_scale[,row.names(subset(ttest_ren_training,p<0.05))]), s = 'lambda.min', type = "response")

dotplot(cbind(Lasso_rep[c(1,3,5,7,9)],Lasso_rep[c(2,4,6,8,10)])[order(cbind(Lasso_rep[c(1,3,5,7,9)],Lasso_rep[c(2,4,6,8,10)])[,1]),], type="b",xlab="Radiomics predictor",ylab="Subject")


ir.pca <- prcomp(rep2021_results_scale, center = FALSE, scale. = FALSE)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1,
              groups = as.factor(c(1,2,1,2,1,2,1,2,1,2)), ellipse = TRUE, var.axes=FALSE,
              circle = TRUE,main="PCA of inter-observer reliability")
g+ geom_point(aes(shape=as.factor(c(1,1,2,2,3,3,4,4,5,5))),size=5,alpha=0.2)
