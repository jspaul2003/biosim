# Title     : biosim
# Objective : investigate distribution of dominant and recessive phenotypes via monte carlo
# Created by: jspaul
# Created on: 2019-04-14
# Required Packages:
#ToDo?: selection at a later date, garbage collector, data representation, debugging

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
    for(a in 1:2){
        for(b in 1:2){
            if(grid[a][b]==1){
                r=r+1
            }
            else if(grid[a][b]==2){
                m=m+1
            }
            else{
                d=d+1
            }
        }
        c(r,m,d)
    }
}
#number of generations,
n=1:1000
#size of population
p=1000
#start simulation of with r only recessive, m recessive and dominant, d dominant. Half of each is female and male approx
r=333
m=333
d=334
#Sexes
r.f=floor(r*0.5)
r.m=r-r.f
m.f=floor(m*0.5)
m.m=m-m.f
d.f=floor(d*0.5)
d.m=d-d.f
#now make the population => EVEN represents male or female, 'tens' digit the type
pop=c()
for(i in 1:r.f){
    pop=c(pop,2)
}
for(i in 1:r.m){
    pop=c(pop,1)
}
for(i in 1:m.f){
    pop=c(pop,12)
}
for(i in 1:m.m){
    pop=c(pop,11)
}
for(i in 1:d.f){
    pop=c(pop,22)
}
for(i in 1:d.m){
    pop=c(pop,21)
}

#begin simulation
for(gen in n){
#loop through the population sort of => tracking reproduction
size=p
chosen=0
    while(chosen<size){
        chosen=chosen+1
        #the chosen one
        ch=floor(runif(1, min=1,max=p+1))
        #is it marked, dont wanna go over it again eg in the 100s
        while(pop[ch]>99){
            ch=floor(runif(1, min=1,max=p+1))
        }
        #the chosen one meets another person-> will they reproduce eg m and f?
        mate=floor(runif(1, min=1,max=p+1))
        #is it marked, dont wanna go over it again eg in the 100s
        while(pop[mate]>99){
            mate=floor(runif(1, min=1,max=p+1))
        }
        if(mate%%2!=ch%%2){
            #success => reproduce
            #we are now less interested in the sexes and more in the tens (recessive, dominant or mix)
            #to do the calculations we will use matrix theory since this is easier
            #p1 parent 1 p2 parent 2
            p1=whoami(pop,ch)
            p2=whoami(pop,mate)
            grid=p1%*%t(p2)
            dis=dist(grid)
            kid=0
            if(runif(1, min=0,max=1)<=dis[0]){
                kid=0
            }
            else if(runif(1, min=0,max=1)>dis[0]+dis[1]){
                kid=20
            }
            else{
                kid=10
            }
            #sex of kid
            if(floor(runif(1, min=1,max=3))==§){
                #male
                kid=kid+1
            }
            else{
            kid=kid+2
            }
            p=+1
            pop=c(pop,kid)
            chosen=chosen+1
            pop[mate]=100+pop[mate]
        }
        pop[ch]=100+pop[ch]
        #and because death exists, theres a 25% chance of death :(
        if(floor(runif(1, min=1,max=5))==1){
            #death
            pop[ch]=0 #represent it as dead
        }
    }
    #new generation, run the garbage collector to reset for next generation, and collect some data!
    newp=0
    newgen =c()
    for(i in 1:p){
        if(pop[p]>100){
            newgen=c(newgen,pop[p]-100)
            newp=newp+1
        }
        else if(pop[p]!=0){
            newgen=c(newgen,pop[p])
            newp=newp+1
        }
    }
    p=newp
    pop=newgen
    #data
}


