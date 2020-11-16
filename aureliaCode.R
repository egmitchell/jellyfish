#import data


##The 2008 network

irishsea2008<-read.delim(file.choose()) #from the inputData2008 tab.
irishsea2008BNI.output<-get.hilo(irishsea2008[,c(1,2,16,19,20,21,22,23)])

#the below file to be run through banjoex 
write.table(irishsea2008BNI.output,file="d:/Dropbox/Projects/y_irishsea/routput/irishsea2008BNI.txt",sep="\t",row.names=FALSE,col.names=FALSE)

#from the command prompt with banjoex and bootstrap installed
banjoex irishsea2008BNI.txt 9 bootstrap=95% repeat=100

#To generate the data for the three decades

cprdataAll<-read.delim(file.choose()) #from the CPRdata70s90s90s tab. includs dates, pci and nao index
temp.all<-read.delim(file.choose()) #from the TempData tab.

species<-as.matrix(cprdataAll[,6:28])
months<-as.matrix(cprdataAll[,5])
years<-as.matrix(cprdataAll[,6])
dates<-as.matrix(cprdataAll[,7]) 
pci<-as.matrix(cprdataAll[,30])
nao<-as.matrix(cprdataAll[,31])

#grouping plankton by biogeographic region cf. Beaugrand 
grouped.data<-matrix(NA,6009,11)
grouped.data[,1]<-species[,22] #warm temperate
grouped.data[,2]<-species[,8] + species[,9] + species[,12] + species[,13] # warm temperate ocean
grouped.data[,3]<-species[,4]+species[,10] + species[,11] #cold temperate
grouped.data[,4]<-species[,1] +species[,2] +species[,3]+species[,20] #shelf sea
grouped.data[,5]<-species[,7]+species[,18] #sub arctic

#working out weekly averages from the input tempdata
req.temp<-cprdataAll[,4:5]

temp.start<-matrix(NA,nrow(temp.all),2)
for(i in 1:nrow(temp.all))
{
	temp.start[i,1]<-(temp.all[i,2]-1971)*12+temp.all[i,1]-1 #gives us month number for 
	temp.start[i,2]<-temp.all[i,3]

}

temp<-matrix(NA,nrow(req.temp),1)
mon<-matrix(NA,nrow(req.temp),1)
for(i in 1:nrow(req.temp))
{
	mon[i,]<-(req.temp[i,2]-1971)*12+req.temp[i,1] #gives us month number for 
	
}

for(i in 1:nrow(req.temp))
{
	temp[i,1]<-mean(temp.start[mon[i,],2])

}

g5<-grouped.data[,1:5]
g5.all<-matrix(NA,nrow(g5),9)
g5.all[,1:5]<-g5
g5.all[,6]<-pci
g5.all[,7]<-nao
g5.all[,8]<-temp
g5.all[,9]<-dates


#need to take weekly averages of everything except years
g5.all.wa<-get.years.for.wa(g5.all,dates) 
g5.all.wa<-remove.zero.row.sp(g5.all.wa,5)

#need to hilo everything except years
g5.hilo.all<-matrix(NA,nrow(g5.all.wa),9)
g5.hilo.all[,1:7]<-get.hilo(g5.all.wa[,1:7])
g5.hilo.all[,8:9]<-get.q(g5.all.wa[,8:9]) #Note that for the NAO and SST we take quartiles not zero high low values 

data.time1<-get.time.period(g5.hilo.all,1971,1980)
data.time2<-get.time.period(g5.hilo.all,1981,1990)
data.time3<-get.time.period(g5.hilo.all,1991,2000)

write.table(data.time1,file="datatime1.txt",sep="\t",row.names=FALSE,col.names=FALSE)
dotify.file(testing.indep(data.time1))

write.table(data.time2,file="datatime2.txt",sep="\t",row.names=FALSE,col.names=FALSE)
dotify.file(testing.indep(data.time2))

write.table(data.time3,file="datatime3.txt",sep="\t",row.names=FALSE,col.names=FALSE)
dotify.file(testing.indep(data.time3))

### Bootstrap and banjoex need to be installed on the computer then the below run. 

banjoex datatime1.txt 10 bootstrap=90% repeat=100
banjoex datatime2.txt 10 bootstrap=90% repeat=100
banjoex datatime3.txt 10 bootstrap=90% repeat=100

#Generating the herring data
envDataHerring<-read.delim(file.choose()) #from the EnvironmentalData&Herring tab.
all.yearly<-get.yearly.averages(g5.all.wa)
herringEnv<-cbind(all.yearly,envDataHerring[,5]])

herringEnv<-cbind(all.yearly,envDataHerring[,5]])

herringEnv.dx<-matrix(NA,nrow(herringEnv),9)
herringEnv.dx[,1:6]<-get.hilo(herringEnv[,1:7])
herringEnv.dx[,1:6]<-get.q(herringEnv[,8:10])#quartiles for sst, nao and herring

write.table(herringEnv.dx[1:9,]file="herring70s.txt",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(herringEnv.dx[10:19,],file="herring80s.txt",sep="\t",row.names=FALSE,col.names=FALSE)
write.table(herringEnv.dx[20:30,],file="herring90s.txt",sep="\t",row.names=FALSE,col.names=FALSE)

#code for working out probablity from: https://zenodo.org/record/3969970#.X7Gfi_P7SUk

##functions used in the above code 
#----------------------------
dotify.file <- function(fp,m)
{
    file.create(paste("c:/",fp,"_avoid.str",sep=""))
    sink(paste("c:/",fp,"avoid.str",sep=""))
    cat(length(m),"\n")
    for (i in 1:length(m))
    {
        cat(m[[i]]);
        cat("\n");
    }
    sink()
}

dotify.indp<-function(filename)
{
	return(dotify.file(filename,testing.indep(filename)))

}



testing.indep<-function(m)
{
    n<-ncol(m)
    m2<-matrix(0,n,n)
    for(i in 1:n)
    {
        for(j in 1:n)
        {
            tbl <- table(m[,i],m[,j])
            m2[i,j] <- ifelse(length(tbl) <= 1 || chisq.test(tbl)$p.value <= 0.25, 1, 0)
        }
    }
   dotify(m2)
   #return(m2)
}


summing.rows<-function(test)
{
	sums<-c(1:nrow(test))
	for(i in 1:nrow(test)) 
	{
		sums[i]<-sum(test[i,])
		
	}
	
	return(sums)
}

summing.rows(test)



matrixify<-function(l)
{
	res<-matrix(NA,length(l),length(l[[1]]))
	for(i in 1:length(l))
	{
		res[i,]<-l[[i]]
	}
	return(res)


}


remove.zero.row<-function(test)
{
	count<-1
	test.list<-list()
	
	if(mean(test[,ncol(test)])<1970)
	{
		sums<-summing.rows(test)
		for(i in 1:nrow(test))
		{
			if(sums[i]!=0)
			{
				test.list[[count]]<-test[i,]
				count<-count+1
			}
		}
	
	}
	if(mean(test[,ncol(test)])>=1970)
	{
		sums<-summing.rows(test[,1:(ncol(test)-1)])
		for(i in 1:nrow(test))
		{
			if(sums[i]!=0)
			{
				test.list[[count]]<-test[i,]
				count<-count+1
			}
		}
	
	}
	
	return(matrixify(test.list))
}



remove.zero.row.sp<-function(data1,sp)#where sp is the number of species ie number of columns to be zeroed
{
	count<-1
	test<-data1[,1:sp]
	test.list<-list()
	
	sums<-summing.rows(test)
	for(i in 1:nrow(test))
	{
		if(sums[i]!=0)
		{
			test.list[[count]]<-data1[i,]
			count<-count+1
		}
	}

	
	
	return(matrixify(test.list))
}


get.median<-function(m)
{
	res<-matrix(NA,1,ncol(m))
	ap<-get.ap(m)
	for(i in 1:ncol(m))
	{
		ap.temp<-as.vector(ap[,i])
		temp<-as.vector(m[,i])
		res[1,i]<-median(temp[ap.temp==1])
	}



	return(res)
}

get.median(data.5g)


get.hilo<-function(m)
{
	meds<-get.median(m)
		
	res<-matrix(NA,nrow(m),ncol(m))
	for(i in 1:ncol(m))
	{
		for(j in 1:nrow(m))
		{
			res[j,i]<-ifelse(m[j,i]==0,0,ifelse(m[j,i]>meds[1,i],2,1))

					}	
	}
	
	return(res)
}

exist.in<-function(m,t1,pm=0)#where end column of m is the years
{
	y<-as.vector(m[,ncol(m)])
	count<-0
	if(pm==1)
	{
		if(t1>max(y))
		{
			return(max(y))
		}
		else
		{
			for(i in 1:length(y))
			{
				if(y[i]==t1)
				{
					count<-count+1
				}
			}
			if(count>0)
			{
				return(t1)
			}
			if(count==0)
			{
				return(exist.in(m,t1+1,1))
			
			}
		}
	}
	if(pm==-1)
	{
		if(t1<min(y))
		{
			return(min(y))
		}
		
		else
		{
			for(i in 1:length(y))
			{
				if(y[i]==t1)
				{
					count<-count+1
				}

			}
			if(count>0)
			{
				return(t1)
			}
			if(count==0)
			{
				return(exist.in(m,t1-1,1))

			}
					
		}
	}
	if(pm==0)
	{
		return("Please input whether to go up 1 or down -1")
	}
}

find.range<-function(m,st,et)
{
	st<-exist.in(m,st,1)
	if(st==0)
	{
		return(c(NA,NA))
	}
	else
	{
		et<-exist.in(m,et,-1)
		return(c(st,et))
	}
	

}


get.time.period<-function(m,start.t,end.t,extra=0)#where final column of m includes years
{
	m<-as.matrix(m)
	r<-find.range(m,start.t,end.t)
	start.t<-r[1]
	end.t<-r[2]
	res<-list()
	count1<-1
	year.list<-m[,ncol(m)]
	for(i in 1:nrow(m))
	{

		if(m[i,ncol(m)]>=start.t)
		{
			if(m[i,ncol(m)]<=end.t)
			{
				res[[count1]]<-m[i,]
				count1<-count1+1
				
			}
		}
	}
	
	res<-matrixify(res)
	if(extra==0)
	{
		res<-res[,1:(ncol(res)-1)]
		return(res)
	}
	if(extra==1)
	{
		res<-res[,1:(ncol(res))]
		return(res)
	}

}

get.yearly.averages<-function(m)#last column should be years
{
	start.year<-m[1,ncol(m)]
	end.year<-m[nrow(m),ncol(m)]
	tot<-end.year-start.year
	count<-0
	res<-matrix(NA,tot,(ncol(m)-1))
	for(i in seq(start.year,(end.year-1),1))
	{
		count<-count+1
		temp<-get.time.period(m,i,i)#
		if(ncol(temp)>(ncol(m)-2))
		{
			temp<-temp[,1:(ncol(temp)-1)]
		}
		
		for(j in 1:ncol(temp))
		{
			#cat(ncol(m)-2,ncol(temp))
			res[count,j]<-mean(temp[,j])
			
		}
		res[count,ncol(temp)+1]<-i
		
	}	
	return(res[1:count,])
}

get.yearly.averages(g5.all.wa)
get.quantile<-function(m)
{
	res<-matrix(NA,nrow(m),ncol(m))
	ap<-get.ap(m)
	for(i in 1:ncol(m))
	{
		ap.temp<-as.vector(ap[,i])
		temp<-as.vector(m[,i])
		res[1:nrow(m),i]<-quantile(temp[ap.temp==1])
	
	
	}



	return(res)
}

get.q<-function(m)
{
	qs<-get.quantile(m)
	
	
	
	res<-matrix(NA,nrow(m),ncol(m))
	for(i in 1:ncol(m))
	{
		#cat(i,"\n")

		for(j in 1:nrow(m))
		{
			#cat("\t",j,"\n")
			res[j,i]<-ifelse(m[j,i]==0,0,ifelse(m[j,i]>qs[4,i],4,ifelse(m[j,i]>qs[3,i],3,ifelse(m[j,i]>qs[2,i],2,1))))


		}	
	}
	
	return(res)
}

get.ap<-function(m)
{
	res<-matrix(NA,nrow(m),ncol(m))
	res<-na.omit(ifelse(m[,1:ncol(m)]>0.00000001, 1, 0)) 

	return(res)
}

get.ap(species)

get.years.for.wa <- function(m,d,extra=0)#d in days from start#extra denotes whether to include month numbers
{
    d<-(as.matrix(d))
    m<-as.matrix(m)
    mres<-matrix(NA,nrow(m),(ncol(m)+1))
    mres[,1:ncol(m)]<-m
    mres[,(ncol(m)+1)]<-floor(d/365)+1971
    m<-mres
    weeks<-get.week.number(d)
    lwn<-list.week.numbers(d)
    
    no.of.weeks<-max(weeks)
    res.no<-count(get.week.number(d))
    res<-list()
    count1<-1
    
    temp<-list()
    for(j in lwn)
    {
        temp<-list()
        count2<-1
        for(i in 1:nrow(d))
        {
            
            if(weeks[i,1]==j)
            {
                temp[[count2]]<-m[i,]
                count2<-count2+1
            }
        }
        temp<-matrixify(temp)
        subres<-matrix(NA,1,ncol(m))
        for(k in 1:ncol(m))
        {
            
            subres[,k]<-mean(temp[,k])
        }
        
        res[[count1]]<-subres
        count1<-count1+1
        
    }
    
    
    return((res))
}

get.week.number <- function(d)
{
    d<-as.matrix(d)
    max.day<-d[nrow(d)]
    no.weeks<-floor(max.day/7)	
    first.day<-d[1,]
    start.row<-1
    week.number<-matrix(NA,nrow(d),1)
    
    for(i in 1:no.weeks)
    {
        start.day<-(i-1)*(7)+first.day
        end.day<-i*7+first.day-1
        #cat(start.day,end.day,"\n")
        for(j in 1:nrow(d))#need to change this to relavant 1:nrow thingy#need to acount incase there are not any data taken with one week
        {
            
            
            if(d[j]>=start.day)
            {
                if(d[j]<=end.day)
                {
                    week.number[j,1]<-i
                    
                }
                
                
                
            }	
        }
        
    }
    
    for(k in 1:nrow(d))
    {
        if(is.na(week.number[k,])==TRUE)
        {
            week.number[k,1]<-no.weeks+1
        }
    }
    
    return(week.number)
}

	