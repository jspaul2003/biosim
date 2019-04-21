# Title     : biosim
# Objective : investigate distribution of dominant and recessive phenotypes via monte carlo
# Created by: jspaul
# Created on: 2019-04-14
# Required Packages:
#ToDo: selection at a later date?, debugging/efficiency?

#NOT SEX Chromosomes
#stuff for later on
whoami = function(pop, p){
    rec=matrix(c(1,1),ncol=1)
    mix=matrix(c(1,2),ncol=1)
    dom=matrix(c(2,2),ncol=1)
    if(pop[p]<5){
        rec
    }
    else if(pop[p]>20){
        dom
    }
    else{
        mix
    }
}
dist = function(grid){
    r=0
    m=0
    d=0
    for(z in 1:4){
        if(grid[z]==1){
            r=r+1
        }
        else if(grid[z]==2){
            m=m+1
        }
        else{
            d=d+1
        }
    }
    c(r,m,d)
}
biosim=function(len,r,m,d,dr){
    data=c()
    a=r
    b=m
    c=d
    #number of generations,
    n=1:len
    #size of population
    p=r+m+d
    #start simulation of with r only recessive, m recessive and dominant, d dominant. Half of each is female and male approx
    #Sexes
    r.f=floor(r*0.5)
    r.m=r-r.f
    m.f=floor(m*0.5)
    m.m=m-m.f
    d.f=floor(d*0.5)
    d.m=d-d.f
    #now make the population => EVEN represents male or female, 'tens' digit the type
    pop=c()
    for(i in 1:(r+m+d)){
        if(r.f!=0){
            pop=c(pop,2)
        }
        if(r.m!=0){
            pop=c(pop,1)
        }
        if(m.f!=0){
            pop=c(pop,12)
        }
        if(m.m!=0){
            pop=c(pop,11)
        }
        if(d.f!=0){
            pop=c(pop,22)
        }
        if(d.f!=0){
            pop=c(pop,21)
        }
    }
    #randomize
    pop=sample(pop)
    print(pop)
    #begin simulation
    for(gen in n){
    #loop through the population sort of => tracking reproduction
    size=p
    chosen=0
        while(chosen<size){
            chosen=chosen+1
            #the chosen one
            ch=floor(runif(1, min=1,max=p+1))
            track=0
            track2=0
            #is it marked, dont wanna go over it again eg in the 100s
            while((pop[ch]>99 || pop[ch]==0) && track<=size){ #what if theyre all dead or done?
                ch=floor(runif(1, min=1,max=p+1))
                track=track+1
                print(track)
                print(size)
                print(pop[ch])
            }
            #the chosen one meets another person-> will they reproduce eg m and f?
            mate=floor(runif(1, min=1,max=p+1))
            #is it marked, dont wanna go over it again eg in the 100s
            while((pop[mate]>99 ||pop[ch]==0) && track2<=size){
                mate=floor(runif(1, min=1,max=p+1))
                track2=track2+1
            }
            #and because death exists, theres a dr chance of death :(
            if(runif(1, min=0,max=1)<=dr){
                #death
                pop[ch]=0 #represent it as dead
            }
            #and because death exists, theres a dr chance of death :(
            if(runif(1, min=0,max=1)<=dr){
                #death
                pop[mate]=0 #represent it as dead
            }
            if(mate%%2!=ch%%2 && pop[mate]!=0 && pop[ch]!=0){
                #success => reproduce
                #we are now less interested in the sexes and more in the tens (recessive, dominant or mix)
                #to do the calculations we will use matrix theory since this is easier
                #p1 parent 1 p2 parent 2
                p1=whoami(pop,ch)
                p2=whoami(pop,mate)
                grid=p1%*%t(p2)
                dis=dist(grid)
                kid=0
                if(runif(1, min=0,max=1)<=dis[1]/4){
                    kid=0
                }
                else if(runif(1, min=0,max=1)<=(dis[3])/4){
                    kid=20

                }
                else{
                    kid=10

                }
                #sex of kid
                if(floor(runif(1, min=1,max=3))==1){
                    #male
                    kid=kid+1
                }
                else{
                kid=kid+2
                }
                if(kid<20){
                    print(pop[ch])
                    print(pop[mate])
                    print(p1)
                    print(p2)
                }
                p=p+1
                pop=c(pop,kid)
                chosen=chosen+1
                pop[mate]=100+pop[mate]
            }
            if(pop[ch]!=0){
                pop[ch]=100+pop[ch]
            }
        }
        #new generation, run the garbage collector to reset for next generation, and collect some data!
        newp=0
        newgen =c()
        if(p<=0){
                break
        }
        else{
            for(i in 1:p){
                if(pop[i]>100){
                    newgen=c(newgen,pop[i]-100)
                    newp=newp+1
                }
                else if(pop[i]!=0){
                    newgen=c(newgen,pop[i])
                    newp=newp+1
                }
            }
        }
        p=newp
        pop=newgen
        pop=sample(pop)
        print(pop)
        r=0
        m=0
        d=0
        if(p<=0){
                break
        }
        else{
            for(i in 1:p){
                if(pop[i]<5){
                    r=r+1
                }
                else if(pop[i]>20){
                    d=d+1
                }
                else{
                    m=m+1
                }
            }
        }
        #data
        a=c(a,r)
        b=c(b,m)
        c=c(c,d)
    }
    data=c(a,b,c)
}
#biosim gets raw data but maybe we want to print that data easily?
easyPrint=function(data){
    par(mfrow=c(2,2))
    plot(data[1:(length(data)/3)],xlab='Generations',ylab='Population of Recessive only carriers')
    plot(data[(1+length(data)/3):(2*length(data)/3)],xlab='Generations',ylab='Population of Mixed carriers')
    plot(data[(1+length(data)/3*2):length(data)],xlab='Generations',ylab='Population of Dominant only carriers')
    recessivenum=c()
    domnum=c()
    for(i in 1:(length(data)/3)){
        print(i)
        print(length(domnum))
        r=data[i]*2+data[i+(length(data)/3)]
        recessivenum=c(recessivenum,r)
        d=data[i+(length(data)/3)]+2*data[i+2*(length(data)/3)]
        domnum=c(domnum,d)
    }
    plot(seq(0,length(data)/3-1),recessivenum,type='l', xlab = 'Generations', ylab='Total Distribution of recessive and dominant genes',col='red')
    legend('bottomright', legend=c("Recessive", "Dominant"), col=c("red", "blue"),pch = c(15,15),bty = "n")
    lines(seq(0,length(data)/3-1),domnum,col='blue')
}


