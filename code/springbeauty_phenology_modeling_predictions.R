# this script is for spring beauty, bee, and canopy closure phenology modeling, future predictions, spring beauty pollination windows calculation and predicted changes


library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(cowplot)


# read spring beauty phenology data
cv.df1 = read_csv(file="Claytonia_virginica_phenodata.csv")

# bee phenology data
ae.df1 = read_csv(file="Andrena_erigeniae_occurence_phenodata.csv")

# canopy phenology data
lsp1 = read_csv(file = "LSP_canopyclosure_data.csv")


# combine data
df_all=rbind.data.frame(cv.df1, ae.df1, lsp1)
df_all$taxa=c(rep("spring beauty", nrow(cv.df1)),
              rep("bee", nrow(ae.df1)),
              rep("canopy", nrow(lsp1)))

# filter ecoregions and check dates
t3=as.data.frame(table(df_all[, c(30,32)]))
t3=t3[t3$Freq>=20,]   # at least 20 records in each ecoregion
t3
t3$NA_L3CODE_taxa=paste0(t3$NA_L3CODE,"_",t3$taxa, sep="")

df1=df_all %>% mutate(NA_L3CODE_taxa=paste0(NA_L3CODE,"_",taxa, sep="")) %>% filter(NA_L3CODE_taxa %in% t3$NA_L3CODE_taxa, NA_L3CODE !="0.0.0")


# model temporal changes -------------------------------------------------------------------------
all.m1 = lmer(doy ~ scale(lat) + scale(elevation)+ scale(year) +
                (1+scale(year)|taxa:NA_L3CODE), 
              data = df1, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(all.m1)
coef(all.m1)

cf1=coef(all.m1)$`taxa:NA_L3CODE`
cf1$NA_L3CODE=sub(".*:", "", rownames(cf1))
cf1$taxa=sub(":.*", "", rownames(cf1))
cf1$taxa_NA_L3CODE=rownames(cf1)
cf1$cf_year=cf1$`scale(year)`/sd(df1$year)    # unstandardized coefficients

cf1 %>% group_by(taxa) %>% summarise(cf_m=mean(cf_year))




# phenological sensitivity modeling-------------------------------------------

# spring beauty ----------------------------
t3=as.data.frame(table(cv.df1$NA_L3CODE))
t3=t3[t3$Freq>=20,]
t3

cv.dfb=cv.df1 %>% filter(NA_L3CODE %in% t3$Var1, NA_L3CODE !="0.0.0") %>% na.omit()   # 21603


# model selection
cv.wave <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tave_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.wave)

cv.wmax <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmax_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.wmax)

cv.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.wmin)

cv.save <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                  (1|NA_L3CODE), 
                data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.save)

cv.smax <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmax_sp_anm)*scale(PPT_sp_anm)+
                  (1|NA_L3CODE), 
                data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.smax)

cv.smin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_sp_anm)*scale(PPT_sp_anm)+
                  (1|NA_L3CODE), 
                data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.smin)

AICc(cv.wave, cv.wmax, cv.wmin, cv.save, cv.smax, cv.smin)

# full model
cv.f <- lmer(doy ~ scale(lat) + scale(elevation)+
               scale(Tmin_wt_anm)*scale(PPT_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
               (1+scale(Tmin_wt_anm)+scale(Tave_sp_anm)|NA_L3CODE), 
             data = cv.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(cv.f)

options(na.action = "na.fail") 
cv.selection = dredge(cv.f)
options(na.action = "na.omit") 

# best model
cv.m <- lmer(doy ~ scale(lat) + scale(elevation)+
               scale(Tmin_wt_anm)+scale(PPT_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
               (1+scale(Tmin_wt_anm)+scale(Tave_sp_anm)|NA_L3CODE), 
             data = cv.dfb, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(cv.m)
coef(cv.m)
r.squaredGLMM(cv.m)

cfcv=coef(cv.m)$NA_L3CODE
cfcv$NA_L3CODE=rownames(cfcv)
cfcv

cv.sds=read_csv("Claytonia_virginica_LT_SD_climate.csv")
colnames(cv.sds)[c(3, 4)]=c("lat", "lon")
cv.sds=cv.sds[, c(-1, -2, -5)]
cv.dfb=left_join(cv.dfb, cv.sds) %>% distinct()

cv.tm=cv.dfb %>% group_by(NA_L3CODE) %>% summarise(mat=mean(LT_MAT, na.rm=T),
                                                   map=mean(LT_MAP, na.rm=T),
                                                   mat_sd=mean(Lsd_MAT, na.rm=T),
                                                   map_sd=mean(Lsd_MAP, na.rm=T),
                                                   tave_isd=mean(isd_Tave, na.rm=T),
                                                   ppt_isd=mean(isd_PPT, na.rm=T),
                                                   
                                                   tm_sp=mean(LT_Tave_sp, na.rm=T),
                                                   pm_sp=mean(LT_PPT_sp, na.rm=T),
                                                   tm_wt=mean(LT_Tmin_wt, na.rm=T),
                                                   pm_wt=mean(LT_PPT_wt, na.rm=T),
                                                   
                                                   tm_sp_sd=mean(Lsd_Tave_sp, na.rm=T),
                                                   pm_sp_sd=mean(Lsd_PPT_sp, na.rm=T),
                                                   tm_sp_isd=mean(isd_Tave_sp, na.rm=T),
                                                   pm_sp_isd=mean(isd_PPT_sp, na.rm=T),
                                                   tm_wt_sd=mean(Lsd_Tmin_wt, na.rm=T),
                                                   pm_wt_sd=mean(Lsd_PPT_wt, na.rm=T),
                                                   tm_wt_isd=mean(isd_Tmin_wt, na.rm=T),
                                                   pm_wt_isd=mean(isd_PPT_wt, na.rm=T))
cv.tm$NA_L3CODE=as.character(cv.tm$NA_L3CODE)
cfcv=left_join(cfcv, cv.tm)
cfcv$cf_us_wt=cfcv$`scale(Tmin_wt_anm)`/sd(cv.dfb$Tmin_wt_anm, na.rm=T)
cfcv$cf_us_sp=cfcv$`scale(Tave_sp_anm)`/sd(cv.dfb$Tave_sp_anm, na.rm=T)


# miner bee --------------------------------
t3=as.data.frame(table(ae.df1$NA_L3CODE))
t3=t3[t3$Freq>=20,]
t3

ae.dfb=ae.df1 %>% filter(NA_L3CODE %in% t3$Var1) 

# model selection
am.wave <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tave_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.wave)

am.wmax <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmax_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.wmax)

am.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.wmin)

am.save <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.save)

am.smax <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmax_sp_anm)*scale(PPT_sp_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.smax)

am.smin <- lmer(doy ~ scale(lat) + scale(elevation)+
                  scale(Tmin_sp_anm)*scale(PPT_sp_anm)+
                  (1|NA_L3CODE), 
                data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.smin)

AICc(am.wave, am.wmax, am.wmin, am.save, am.smax, am.smin)


# full model
am.f <- lmer(doy ~ scale(lat) + scale(elevation)+
               scale(Tmin_wt_anm)*scale(PPT_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
               (1+scale(Tave_sp_anm)|NA_L3CODE), 
             data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.f) 

options(na.action = "na.fail") 
am.selection=dredge(am.f)
options(na.action = "na.omit") 

# the best model
am.m1 <- lmer(doy ~ scale(lat) + scale(elevation)+
                scale(Tmin_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                (1+scale(Tave_sp_anm)|NA_L3CODE), 
              data = ae.dfb, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(am.m1)
coef(am.m1)
r.squaredGLMM(am.m1)

cfam1=coef(am.m1)$NA_L3CODE
cfam1$NA_L3CODE=rownames(cfam1)
cfam1

ae.sds=read_csv("Andrena_erigeniae_LT_SD_climate.csv")
colnames(ae.sds)[c(3, 4)]=c("lat", "lon")
ae.sds=ae.sds[, c(-1, -2, -5)]
ae.dfb=left_join(ae.dfb, ae.sds)

am.tm=ae.dfb %>% group_by(NA_L3CODE) %>% summarise(mat=mean(LT_MAT, na.rm=T),
                                                   map=mean(LT_MAP, na.rm=T),
                                                   mat_sd=mean(Lsd_MAT, na.rm=T),
                                                   map_sd=mean(Lsd_MAP, na.rm=T),
                                                   tave_isd=mean(isd_Tave, na.rm=T),
                                                   ppt_isd=mean(isd_PPT, na.rm=T),
                                                   
                                                   tm_sp=mean(LT_Tave_sp, na.rm=T),
                                                   pm_sp=mean(LT_PPT_sp, na.rm=T),
                                                   tm_wt=mean(LT_Tmin_wt, na.rm=T),
                                                   pm_wt=mean(LT_PPT_wt, na.rm=T),
                                                   
                                                   tm_sp_sd=mean(Lsd_Tave_sp, na.rm=T),
                                                   pm_sp_sd=mean(Lsd_PPT_sp, na.rm=T),
                                                   tm_sp_isd=mean(isd_Tave_sp, na.rm=T),
                                                   pm_sp_isd=mean(isd_PPT_sp, na.rm=T),
                                                   tm_wt_sd=mean(Lsd_Tmin_wt, na.rm=T),
                                                   pm_wt_sd=mean(Lsd_PPT_wt, na.rm=T),
                                                   tm_wt_isd=mean(isd_Tmin_wt, na.rm=T),
                                                   pm_wt_isd=mean(isd_PPT_wt, na.rm=T))
am.tm$NA_L3CODE=as.character(am.tm$NA_L3CODE)
cfam1=left_join(cfam1, am.tm)
cfam1$cf_us_sp=cfam1$`scale(Tave_sp_anm)`/sd(ae.dfb$Tave_sp_anm, na.rm=T)


# full model
am.f <- lmer(doy ~ scale(lat) + scale(elevation)+
               scale(Tmin_wt_anm)*scale(PPT_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
               (1+scale(Tmin_wt_anm)|NA_L3CODE), 
             data = ae.dfb, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(am.f) 

options(na.action = "na.fail") 
am.selection=dredge(am.f)
options(na.action = "na.omit") 

# the best model
am.m2 <- lmer(doy ~ scale(lat) + scale(elevation)+
                scale(Tmin_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                (1+scale(Tmin_wt_anm)||NA_L3CODE), 
              data = ae.dfb, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(am.m2)
coef(am.m2)
r.squaredGLMM(am.m2)

cfam2=coef(am.m2)$NA_L3CODE
cfam2$NA_L3CODE=rownames(cfam2)
cfam2

cfam2=left_join(cfam2, am.tm)
cfam2$cf_us_wt=cfam2$`scale(Tmin_wt_anm)`/sd(ae.dfb$Tmin_wt_anm, na.rm=T)

ggplot(cfam2, aes(x=tm_wt, y=cf_us_wt))+
  geom_point()+
  geom_smooth(method="lm", se=F)

cor.test(cfam2$tm_wt, cfam2$cf_us_wt)  


# canopy closure -------------------------
t3=as.data.frame(table(lsp1$NA_L3CODE))
t3=t3[t3$Freq>=20,]
t3

lsp1a=lsp1 %>% filter(NA_L3CODE %in% t3$Var1, NA_L3CODE !="0.0.0")   # 315749


# model selection
lsp.wave <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tave_wt_anm)*scale(PPT_wt_anm)+
                   (1|NA_L3CODE), 
                 data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.wave)

lsp.wmax <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmax_wt_anm)*scale(PPT_wt_anm)+
                   (1|NA_L3CODE), 
                 data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.wmax)

lsp.wmin <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmin_wt_anm)*scale(PPT_wt_anm)+
                   (1|NA_L3CODE), 
                 data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.wmin)

lsp.save <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                   (1|NA_L3CODE), 
                 data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.save)

lsp.smax <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmax_sp_anm)*scale(PPT_sp_anm)+
                   (1|NA_L3CODE), 
                 data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.smax)

lsp.smin <- lmer(doy ~ scale(lat) + scale(elevation)+
                   scale(Tmin_sp_anm)*scale(PPT_sp_anm)+
                   (1|NA_L3CODE), 
                 data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.smin)

AICc(lsp.wave, lsp.wmax, lsp.wmin, lsp.save, lsp.smax, lsp.smin)

# full model
lsp.f <- lmer(doy ~ scale(lat) + scale(elevation)+
                scale(Tmin_wt_anm)*scale(PPT_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                (1+scale(Tmin_wt_anm)+scale(Tave_sp_anm)|NA_L3CODE), 
              data = lsp1a, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(lsp.f)

options(na.action = "na.fail") 
lsp.selection = dredge(lsp.f)
options(na.action = "na.omit") 

# best model
lsp.m <- lmer(doy ~ scale(lat) + scale(elevation)+
                scale(Tmin_wt_anm)*scale(PPT_wt_anm)+scale(Tave_sp_anm)*scale(PPT_sp_anm)+
                (1+scale(Tmin_wt_anm)+scale(Tave_sp_anm)|NA_L3CODE), 
              data = lsp1a, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(lsp.m)
coef(lsp.m)
r.squaredGLMM(lsp.m)


cflsp=coef(lsp.m)$NA_L3CODE
cflsp$NA_L3CODE=rownames(cflsp)
cflsp

lsp.sds = rbind.data.frame(cv.sds, ae.sds) %>% distinct()
lsp1a=left_join(lsp1a, lsp.sds)

cv.tm=lsp1a %>% group_by(NA_L3CODE) %>% summarise(mat=mean(LT_MAT, na.rm=T),
                                                  map=mean(LT_MAP, na.rm=T),
                                                  mat_sd=mean(Lsd_MAT, na.rm=T),
                                                  map_sd=mean(Lsd_MAP, na.rm=T),
                                                  tave_isd=mean(isd_Tave, na.rm=T),
                                                  ppt_isd=mean(isd_PPT, na.rm=T),
                                                  
                                                  tm_sp=mean(LT_Tave_sp, na.rm=T),
                                                  pm_sp=mean(LT_PPT_sp, na.rm=T),
                                                  tm_wt=mean(LT_Tmin_wt, na.rm=T),
                                                  pm_wt=mean(LT_PPT_wt, na.rm=T),
                                                  
                                                  tm_sp_sd=mean(Lsd_Tave_sp, na.rm=T),
                                                  pm_sp_sd=mean(Lsd_PPT_sp, na.rm=T),
                                                  tm_sp_isd=mean(isd_Tave_sp, na.rm=T),
                                                  pm_sp_isd=mean(isd_PPT_sp, na.rm=T),
                                                  tm_wt_sd=mean(Lsd_Tmin_wt, na.rm=T),
                                                  pm_wt_sd=mean(Lsd_PPT_wt, na.rm=T),
                                                  tm_wt_isd=mean(isd_Tmin_wt, na.rm=T),
                                                  pm_wt_isd=mean(isd_PPT_wt, na.rm=T))
cv.tm$NA_L3CODE=as.character(cv.tm$NA_L3CODE)
cflsp=left_join(cflsp, cv.tm)
cflsp$cf_us_wt=cflsp$`scale(Tmin_wt_anm)`/sd(lsp1a$Tmin_wt_anm)
cflsp$cf_us_sp=cflsp$`scale(Tave_sp_anm)`/sd(lsp1a$Tave_sp_anm)


# model predictions ----------------------------------------------------------
# current year estimation (1981-2010)

# calculate 30y normal ecoregion level climate variables
cv.sy=read_csv("springbeauty_locations_all_by2022_20240103_1901-2022SY.csv")
ae.sy=read_csv("minerbee_locations_all_20240103_1901-2022SY.csv") 

# cv.loc=read_csv(file="data_raw/locations/springbeauty_locations_all_20231219.csv")
# summary(cv.loc)

cv.LT= cv.sy %>% rename(lat=Latitude, lon=Longitude, elevation=Elevation) %>% 
  group_by(lat, lon, elevation) %>% 
  summarise(LT_Tmin_wt=mean(Tmin_wt), LT_Tave_sp=mean(Tave_sp),
            LT_PPT_wt=mean(PPT_wt),
            LT_PPT_sp=mean(PPT_sp)) %>% ungroup()

colnames(cv.sy)

cv.sy1 = cv.sy[, c(1, 4:6, 11, 16, 19, 20)] %>% 
  rename(year=Year, lat=Latitude, lon=Longitude, elevation=Elevation) %>%
  filter(year>=1981, year<=2010) 



ae.LT= ae.sy %>% rename(lat=Latitude, lon=Longitude, elevation=Elevation) %>% 
  group_by(lat, lon, elevation) %>% 
  summarise(LT_Tmin_wt=mean(Tmin_wt), LT_Tave_sp=mean(Tave_sp),
            LT_PPT_wt=mean(PPT_wt),
            LT_PPT_sp=mean(PPT_sp)) %>% ungroup()

ae.sy1 = ae.sy[, c(1, 4:6, 11, 16, 19, 20)] %>% 
  rename(year=Year, lat=Latitude, lon=Longitude, elevation=Elevation) %>%
  filter(year>=1981, year<=2010) 


all.LT= rbind.data.frame(ungroup(cv.LT), ungroup(ae.LT)) %>% distinct()



# prepare new data
# spatial and random effect groups
cv.sp=cv.dfb %>% select(lat, lon, elevation, NA_L3CODE) %>% distinct()
ae.sp=ae.dfb %>% select(lat, lon, elevation, NA_L3CODE) %>% distinct()
#lsp.sp=rbind.data.frame(cv.sp, ae.sp) %>% distinct()

# climate
cv.clm = cv.sy1 %>% group_by(lat, lon) %>% summarise(Tmwt=mean(Tmin_wt), Tmsp=mean(Tave_sp), 
                                                     Pmwt=mean(PPT_wt), Pmsp=mean(PPT_sp))
cv.clm = left_join(cv.clm, cv.LT) %>% distinct()

cv.clm = cv.clm %>% mutate(Tmin_wt_anm = Tmwt-LT_Tmin_wt, 
                           Tave_sp_anm = Tmsp-LT_Tave_sp,
                           PPT_wt_anm = (Pmwt-LT_PPT_wt)/LT_PPT_wt,
                           PPT_sp_anm = (Pmsp-LT_PPT_sp)/LT_PPT_sp)

cv.nd = left_join(cv.sp, cv.clm[, c(1:2, 11:14)])


ae.clm = ae.sy1 %>% group_by(lat, lon) %>% summarise(Tmwt=mean(Tmin_wt), Tmsp=mean(Tave_sp), 
                                                     Pmwt=mean(PPT_wt), Pmsp=mean(PPT_sp))
ae.clm = left_join(ae.clm, ae.LT) %>% distinct()

ae.clm = ae.clm %>% mutate(Tmin_wt_anm = Tmwt-LT_Tmin_wt, 
                           Tave_sp_anm = Tmsp-LT_Tave_sp,
                           PPT_wt_anm = (Pmwt-LT_PPT_wt)/LT_PPT_wt,
                           PPT_sp_anm = (Pmsp-LT_PPT_sp)/LT_PPT_sp)

ae.nd = left_join(ae.sp, ae.clm[, c(1:2, 11:14)])

lsp.nd = rbind.data.frame(ae.nd, cv.nd) %>% distinct()

all.nd = rbind.data.frame(ae.nd, cv.nd) %>% distinct() %>% filter(NA_L3CODE %in% ae.dfb$NA_L3CODE)

# 30-y average
all.nd$curr.cv=predict(cv.m, newdata=all.nd)
all.nd$curr.ae1=predict(am.m1, newdata=all.nd)
all.nd$curr.ae2=predict(am.m2, newdata=all.nd)
all.nd$curr.lsp=predict(lsp.m, newdata=all.nd)

# future prediction
cv.gcm4=read_csv("data_raw/locations/springbeauty_locations_all_by2022_20240103_4 GCMsSY.csv") 

cv.gcm4= cv.gcm4 %>% 
  rename(year=Year, lat=Latitude, lon=Longitude, elevation=Elevation) %>% 
  select(year, lat , lon, elevation, Tmin_wt, Tave_sp, PPT_wt, PPT_sp) 

unique(cv.gcm4$year)

cv.gcm4=left_join(cv.gcm4, cv.LT) %>% 
  mutate(Tmin_wt_anm = Tmin_wt-LT_Tmin_wt, 
         Tave_sp_anm = Tave_sp-LT_Tave_sp, 
         PPT_wt_anm = (PPT_wt-LT_PPT_wt)/LT_PPT_wt,
         PPT_sp_anm = (PPT_sp-LT_PPT_sp)/LT_PPT_sp)


# 13GCM-ensemble-ssp245-2041-2060
cv.fpd1 = left_join(cv.nd[, c(1:4)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp245_2041-2060.gcm")) %>% distinct()
cv.nd$ssp245_4160=predict(cv.m, newdata=cv.fpd1)

# 13GCM-ensemble-ssp585-2041-2060
cv.fpd2 = left_join(cv.nd[, c(1:4)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp585_2041-2060.gcm")) %>% distinct()
cv.nd$ssp585_4160=predict(cv.m, newdata=cv.fpd2)

# 13GCM-ensemble-ssp245-2081-2100
cv.fpd3 = left_join(cv.nd[, c(1:4)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp245_2081-2100.gcm")) %>% distinct()
cv.nd$ssp245_8100=predict(cv.m, newdata=cv.fpd3)

# 13GCM-ensemble-ssp585-2081-2100
cv.fpd4 = left_join(cv.nd[, c(1:4)], filter(cv.gcm4, year=="13GCMs_ensemble_ssp585_2081-2100.gcm")) %>% distinct()
cv.nd$ssp585_8100=predict(cv.m, newdata=cv.fpd4)

summary(cv.nd)


# spring beauty miner bee
ae.gcm4=read_csv("data_raw/locations/minerbee_locations_all_20240103_4 GCMsSY.csv") 

ae.gcm4 = ae.gcm4 %>% 
  rename(year=Year, lat=Latitude, lon=Longitude, elevation=Elevation) %>% 
  dplyr::select(year, lat , lon, elevation, Tmin_wt, Tave_sp, PPT_wt, PPT_sp) 


unique(ae.gcm4$year)

ae.gcm4=left_join(ae.gcm4, ae.LT) %>% 
  mutate(Tmin_wt_anm = Tmin_wt-LT_Tmin_wt, 
         Tave_sp_anm = Tave_sp-LT_Tave_sp, 
         PPT_wt_anm = (PPT_wt-LT_PPT_wt)/LT_PPT_wt,
         PPT_sp_anm = (PPT_sp-LT_PPT_sp)/LT_PPT_sp)


ae.gcm41 = left_join(ae.gcm4, distinct(ae.dfb[, c(1,2,6,52)]))

ae.gcm41 %>% na.omit() %>% group_by(year, NA_L3CODE) %>% 
  summarise(Tw_m=mean(Tmin_wt_anm, na.rm=T), Ts_m=mean(Tave_sp_anm, na.rm=T),
            Pw_m=mean(PPT_wt_anm, na.rm=T), Ps_m=mean(PPT_sp_anm, na.rm=T)) %>% view()


# 13GCM-ensemble-ssp245-2041-2060
ae.fpd1 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp245_2041-2060.gcm")) %>% distinct()
ae.nd$ssp245_4160_1=predict(am.m1, newdata=ae.fpd1)
ae.nd$ssp245_4160_2=predict(am.m2, newdata=ae.fpd1)

# 13GCM-ensemble-ssp585-2041-2060
ae.fpd2 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp585_2041-2060.gcm")) %>% distinct()
ae.nd$ssp585_4160_1=predict(am.m1, newdata=ae.fpd2)
ae.nd$ssp585_4160_2=predict(am.m2, newdata=ae.fpd2)

# 13GCM-ensemble-ssp245-2081-2100
ae.fpd3 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp245_2081-2100.gcm")) %>% distinct()
ae.nd$ssp245_8100_1=predict(am.m1, newdata=ae.fpd3)
ae.nd$ssp245_8100_2=predict(am.m2, newdata=ae.fpd3)

# 13GCM-ensemble-ssp585-2081-2100
ae.fpd4 = left_join(ae.nd[, c(1:4)], filter(ae.gcm4, year=="13GCMs_ensemble_ssp585_2081-2100.gcm")) %>% distinct()
ae.nd$ssp585_8100_1=predict(am.m1, newdata=ae.fpd4)
ae.nd$ssp585_8100_2=predict(am.m2, newdata=ae.fpd4)

summary(ae.nd)


# 13GCM-ensemble-ssp245-2041-2060
lsp.fpd1 = rbind.data.frame(cv.fpd1, ae.fpd1) %>% distinct()
lsp.nd$ssp245_4160=predict(lsp.m, newdata=lsp.fpd1)

# 13GCM-ensemble-ssp585-2041-2060
lsp.fpd2 = rbind.data.frame(cv.fpd2, ae.fpd2) %>% distinct()
lsp.nd$ssp585_4160=predict(lsp.m, newdata=lsp.fpd2)

# 13GCM-ensemble-ssp245-2081-2100
lsp.fpd3 = rbind.data.frame(cv.fpd3, ae.fpd3) %>% distinct()
lsp.nd$ssp245_8100=predict(lsp.m, newdata=lsp.fpd3)

# 13GCM-ensemble-ssp585-2081-2100
lsp.fpd4 = rbind.data.frame(cv.fpd4, ae.fpd4) %>% distinct()
lsp.nd$ssp585_8100=predict(lsp.m, newdata=lsp.fpd4)

summary(lsp.nd)


all.fpd1 = rbind.data.frame(cv.fpd1, ae.fpd1) %>% distinct() %>% filter(NA_L3CODE %in% ae.dfb$NA_L3CODE)
all.nd$ssp245_4160.cv=predict(cv.m, newdata=all.fpd1)
all.nd$ssp245_4160.ae1=predict(am.m1, newdata=all.fpd1)
all.nd$ssp245_4160.ae2=predict(am.m2, newdata=all.fpd1)
all.nd$ssp245_4160.lsp=predict(lsp.m, newdata=all.fpd1)


all.fpd2 = rbind.data.frame(cv.fpd2, ae.fpd2) %>% distinct() %>% filter(NA_L3CODE %in% ae.dfb$NA_L3CODE)
all.nd$ssp585_4160.cv=predict(cv.m, newdata=all.fpd2)
all.nd$ssp585_4160.ae1=predict(am.m1, newdata=all.fpd2)
all.nd$ssp585_4160.ae2=predict(am.m2, newdata=all.fpd2)
all.nd$ssp585_4160.lsp=predict(lsp.m, newdata=all.fpd2)


all.fpd3 = rbind.data.frame(cv.fpd3, ae.fpd3) %>% distinct() %>% filter(NA_L3CODE %in% ae.dfb$NA_L3CODE)
all.nd$ssp245_8100.cv=predict(cv.m, newdata=all.fpd3)
all.nd$ssp245_8100.ae1=predict(am.m1, newdata=all.fpd3)
all.nd$ssp245_8100.ae2=predict(am.m2, newdata=all.fpd3)
all.nd$ssp245_8100.lsp=predict(lsp.m, newdata=all.fpd3)


all.fpd4 = rbind.data.frame(cv.fpd4, ae.fpd4) %>% distinct() %>% filter(NA_L3CODE %in% ae.dfb$NA_L3CODE)
all.nd$ssp585_8100.cv=predict(cv.m, newdata=all.fpd4)
all.nd$ssp585_8100.ae1=predict(am.m1, newdata=all.fpd4)
all.nd$ssp585_8100.ae2=predict(am.m2, newdata=all.fpd4)
all.nd$ssp585_8100.lsp=predict(lsp.m, newdata=all.fpd4)


# same location predictions all.nd
scn5=c("fit.curr","ssp245_4160","ssp245_8100","ssp585_4160","ssp585_8100")

cv.pd=all.nd[, c(1:4, 9, 13, 17, 21, 25)]
colnames(cv.pd)[5:9]=scn5
ae.pd1=all.nd[, c(1:4, 10, 14, 18, 22, 26)]
colnames(ae.pd1)[5:9]=scn5
ae.pd2=all.nd[, c(1:4, 11, 15, 19, 23, 27)]
colnames(ae.pd2)[5:9]=scn5
lsp.pd=all.nd[, c(1:4, 12, 16, 20, 24, 28)]
colnames(lsp.pd)[5:9]=scn5




# calculate pollination period (end-start)
# end - the latest date, usually the canopy closure date
# start - the latter date between spring beauty and bee
colnames(all.nd)

pl_curr = all.nd %>% dplyr::select(1:4, curr.cv:curr.lsp)

pl_curr$start1 = pl_curr$curr.ae1
pl_curr$start1[which(pl_curr$curr.ae1<pl_curr$curr.cv)]=pl_curr$curr.cv[which(pl_curr$curr.ae1<pl_curr$curr.cv)]

pl_curr$end = pl_curr$curr.lsp
pl_curr$end[pl_curr$end<pl_curr$curr.ae1]  

pl_curr$pollination1=pl_curr$end - pl_curr$start1
pl_curr$peak=as.integer((pl_curr$start1+pl_curr$end)/2)

pl_curr %>% group_by(NA_L3CODE) %>% summarise(pl.m=mean(pollination1, na.rm=T))

ggplot(pl_curr, aes(x=NA_L3CODE, y=pollination1))+
  geom_boxplot()


pl_ssp245a= all.nd %>% dplyr::select(1:4, ssp245_4160.cv:ssp245_4160.lsp)

pl_ssp245a$start1 = pl_ssp245a$ssp245_4160.ae1
pl_ssp245a$start1[which(pl_ssp245a$ssp245_4160.ae1<pl_ssp245a$ssp245_4160.cv)]=pl_ssp245a$ssp245_4160.cv[which(pl_ssp245a$ssp245_4160.ae1<pl_ssp245a$ssp245_4160.cv)]

pl_ssp245a$end = pl_ssp245a$ssp245_4160.lsp
pl_ssp245a$end[pl_ssp245a$end<pl_ssp245a$ssp245_4160.ae1]  

pl_ssp245a$pollination1=pl_ssp245a$end - pl_ssp245a$start1
pl_ssp245a$peak=as.integer((pl_ssp245a$start1+pl_ssp245a$end)/2)

pl_ssp245a %>% group_by(NA_L3CODE) %>% summarise(pl.m=mean(pollination1, na.rm=T))

ggplot(pl_ssp245a, aes(x=NA_L3CODE, y=pollination1))+
  geom_boxplot()


pl_ssp245b= all.nd %>% dplyr::select(1:4, ssp245_8100.cv:ssp245_8100.lsp)

pl_ssp245b$start1 = pl_ssp245b$ssp245_8100.ae1
pl_ssp245b$start1[which(pl_ssp245b$ssp245_8100.ae1<pl_ssp245b$ssp245_8100.cv)]=pl_ssp245b$ssp245_8100.cv[which(pl_ssp245b$ssp245_8100.ae1<pl_ssp245b$ssp245_8100.cv)]

pl_ssp245b$end = pl_ssp245b$ssp245_8100.lsp
pl_ssp245b$end[pl_ssp245b$end<pl_ssp245b$ssp245_8100.ae1]  

pl_ssp245b$pollination1=pl_ssp245b$end - pl_ssp245b$start1
pl_ssp245b$peak=as.integer((pl_ssp245b$start1+pl_ssp245b$end)/2)


pl_ssp585a= all.nd %>% dplyr::select(1:4, ssp585_4160.cv:ssp585_4160.lsp)

pl_ssp585a$start1 = pl_ssp585a$ssp585_4160.ae1
pl_ssp585a$start1[which(pl_ssp585a$ssp585_4160.ae1<pl_ssp585a$ssp585_4160.cv)]=pl_ssp585a$ssp585_4160.cv[which(pl_ssp585a$ssp585_4160.ae1<pl_ssp585a$ssp585_4160.cv)]

pl_ssp585a$end = pl_ssp585a$ssp585_4160.lsp
pl_ssp585a$end[pl_ssp585a$end<pl_ssp585a$ssp585_4160.ae1]  

pl_ssp585a$pollination1=pl_ssp585a$end - pl_ssp585a$start1
pl_ssp585a$peak=as.integer((pl_ssp585a$start1+pl_ssp585a$end)/2)

pl_ssp585a %>% group_by(NA_L3CODE) %>% summarise(pl.m=mean(pollination1, na.rm=T))

ggplot(pl_ssp585a, aes(x=NA_L3CODE, y=pollination1))+
  geom_boxplot()


pl_ssp585b= all.nd %>% dplyr::select(1:4, ssp585_8100.cv:ssp585_8100.lsp)

pl_ssp585b$start1 = pl_ssp585b$ssp585_8100.ae1
pl_ssp585b$start1[which(pl_ssp585b$ssp585_8100.ae1<pl_ssp585b$ssp585_8100.cv)]=pl_ssp585b$ssp585_8100.cv[which(pl_ssp585b$ssp585_8100.ae1<pl_ssp585b$ssp585_8100.cv)]

pl_ssp585b$end = pl_ssp585b$ssp585_8100.lsp
pl_ssp585b$end[pl_ssp585b$end<pl_ssp585b$ssp585_8100.ae1]  

pl_ssp585b$pollination1=pl_ssp585b$end - pl_ssp585b$start1
pl_ssp585b$peak=as.integer((pl_ssp585b$start1+pl_ssp585b$end)/2)


colnames(pl_curr)[5:8]=c("cv", "ae1", "ae2", "lsp")
colnames(pl_ssp245a)[5:8]=c("cv", "ae1", "ae2", "lsp")
colnames(pl_ssp245b)[5:8]=c("cv", "ae1", "ae2", "lsp")
colnames(pl_ssp585a)[5:8]=c("cv", "ae1", "ae2", "lsp")
colnames(pl_ssp585b)[5:8]=c("cv", "ae1", "ae2", "lsp")

pl_all=rbind.data.frame(pl_curr, pl_ssp245a, pl_ssp245b, pl_ssp585a, pl_ssp585b)
pl_all$scenario = rep(c("curr","ssp245_4160","ssp245_8100","ssp585_4160","ssp585_8100"),each=nrow(pl_curr))


# changes in pollination period and mean date

pl_ssp245a$d_period=pl_ssp245a$pollination1 - pl_curr$pollination1
pl_ssp245a$d_date=pl_ssp245a$peak - pl_curr$peak
pl_ssp245a$p_period=(pl_ssp245a$pollination1 - pl_curr$pollination1)/pl_curr$pollination1
pl_ssp245a$p_date=(pl_ssp245a$peak - pl_curr$peak)/pl_curr$peak


pl_ssp245b$d_period=pl_ssp245b$pollination1 - pl_curr$pollination1
pl_ssp245b$d_date=pl_ssp245b$peak - pl_curr$peak
pl_ssp245b$p_period=(pl_ssp245b$pollination1 - pl_curr$pollination1)/pl_curr$pollination1
pl_ssp245b$p_date=(pl_ssp245b$peak - pl_curr$peak)/pl_curr$peak


pl_ssp585a$d_period=pl_ssp585a$pollination1 - pl_curr$pollination1
pl_ssp585a$d_date=pl_ssp585a$peak - pl_curr$peak
pl_ssp585a$p_period=(pl_ssp585a$pollination1 - pl_curr$pollination1)/pl_curr$pollination1
pl_ssp585a$p_date=(pl_ssp585a$peak - pl_curr$peak)/pl_curr$peak


pl_ssp585b$d_period=pl_ssp585b$pollination1 - pl_curr$pollination1
pl_ssp585b$d_date=pl_ssp585b$peak - pl_curr$peak
pl_ssp585b$p_period=(pl_ssp585b$pollination1 - pl_curr$pollination1)/pl_curr$pollination1
pl_ssp585b$p_date=(pl_ssp585b$peak - pl_curr$peak)/pl_curr$peak


pl_ch=rbind.data.frame(pl_ssp245a, pl_ssp245b, pl_ssp585a, pl_ssp585b)
pl_ch$scenario = rep(c("ssp245_4160","ssp245_8100","ssp585_4160","ssp585_8100"),each=nrow(pl_ssp245a))
