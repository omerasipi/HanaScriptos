import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class LebesgueIntegrate:
    def __init__(self, f, a, b, n, m=400, fmin=None, fmax=None):
        import warnings
        warnings.filterwarnings("ignore")

        self.n = n
        self.m = np.max([m,1000*n])
        self.xp = np.linspace(a, b, self.m)
        self.yp = f(self.xp)
        
        ind = np.isfinite(self.yp)
        self.yp = self.yp[ind]
        self.xp = self.xp[ind]
        self.m = self.yp.shape[0]

        if fmin==None:
            fmin = np.min(self.yp)-1e-14
        if fmax==None:
            fmax = np.max(self.yp)+1e-14
        self.yi = np.linspace(fmin,fmax,self.n+1)
        
    def computeSupport(self, k):
        ind = []
        ind2 = []
        for i in range(self.m):
            if (self.yi[k] < self.yp[i]) and (self.yp[i] <= self.yi[k+1]):
                ind.append(i)
                if len(ind)>1:
                    if self.xp[i]>self.xp[ind[-2]+1]:
                        ind2.append(ind[-2])
                        ind2.append(i)
                else:
                    ind2.append(i)
        if len(ind)>0:
            ind2.append(ind[-1])
            ind2=np.array(ind2).reshape((int(len(ind2)/2),2))
        return ind, ind2
    
    def computeMeasure(self,ind2):
        meas = 0.
        if len(ind2)>0:
            for i,j in ind2:
                meas += self.xp[j]-self.xp[i]
        return meas

    def integrate(self):
        Lmin = 0
        Lmax = 0
        for k in range(self.n):
            ind, ind2 = self.computeSupport(k)
            meas = self.computeMeasure(ind2)
            Lmin += meas*self.yi[k]
            Lmax += meas*self.yi[k+1]
        return Lmin, Lmax

    def myPlot(self, sel, ylim=None):
        colors=list(mcolors.TABLEAU_COLORS.keys())
        if len(colors) < self.n:
            colors *= int(np.ceil(self.n/len(colors)))
        if ylim == None:
            ylim = (self.yi[0]-0.05*(self.yi[-1]-self.yi[0]),self.yi[-1]+0.05*(self.yi[-1]-self.yi[0]))
        plt.ylim(ylim)
        plt.plot(self.xp, self.yp,'-')
        for k,col in zip(range(self.n),colors):
            ind, ind2 = self.computeSupport(k)
            for i,j in ind2:
                plt.plot(self.xp[range(i,j)],self.yp[range(i,j)],c=col)
                plt.plot(self.xp[range(i,j)],(self.yi[0]-0.025*(self.yi[-1]-self.yi[0]))*np.ones(j-i),c=col)#,lw=4)
        ind, ind2 = self.computeSupport(sel)
        for i,j in ind2:
            plt.fill_between(self.xp[range(i,j)],self.yp[range(i,j)],color=colors[sel],alpha=0.3)
        plt.fill_betweenx(self.yi[[sel,sel+1]],self.xp[[-1,-1]],color=colors[sel],alpha=0.15)
        for i in ind2.flatten():
            plt.axvline(self.xp[i], c='tab:gray',lw=0.5)
        for yii in self.yi:
            plt.axhline(yii, c='tab:gray',lw=0.5)
        plt.title(r'Illustration Lebesgue Integral $\mu(E_'+str(sel)+') = '+str(np.round(self.computeMeasure(ind2),3))+'$')
        plt.show()

