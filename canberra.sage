class singular_problem:
    def __init__(self, x, f, g):
        self.x = SR(x)
        self.f = SR(f)
        self.g = SR(g)
    def latex(self):
        print("\\left \\{ \\begin{aligned} &")
        print("".join([r'\mu\frac{d'+latex(self.x)+'}{dt}'+'='+latex(self.f) +', \\\\ &']))
        print("".join([latex(self.x)+'(0)='+latex(self.g)]))
        print("\\end{aligned} \\right. ")
    def r_series(self,expr,n):
        if expr==self.x:
            return sum([function('r'+str(expr)+str(nn))(t)*mu^nn for nn in range(n+1)])
        else:
            pr=singular_problem(self.x,self.f,self.g)
            return SR(expr).subs(self.x==pr.r_series(self.x,n)).series(mu,n+1).truncate()
    def r_list(self,n):
        pr=singular_problem(self.x,self.f,self.g)
        return pr.r_series(self.x,n).coefficients(sparse=False)
    def r_odes(self,n):
        pr=singular_problem(self.x,self.f,self.g)
        ode=mu*diff(pr.r_series(self.x, n-1),t)-pr.r_series(self.f, n)
        return ode.coefficients(sparse=False)
    def r_solve(self,n):
        pr=singular_problem(self.x,self.f,self.g)
        eqs=pr.r_odes(n)
        R=pr.r_list(n)
        sol=[]
        k=0
        for (eq,r) in zip(eqs,R):
            eq=eq.subs(sol)
            for nn in range(n+1):
                if k==0:
                    eq=diff(eq,t).subs(sol)
                    k=1
                else:
                    sol.append(diff(r,t,nn)==-eq.subs(diff(r,t,nn)==0)/eq.coefficient(diff(r,t,nn)))
                    if nn<n:
                        eq=diff(eq,t).subs(sol)
        return sol
    
    def p_series(self,expr,n):
        if expr==self.x:
            return sum([function('p'+str(expr)+str(nn))(tau)*mu^nn for nn in range(n+1)])
        else:
            pr=singular_problem(self.x,self.f,self.g)
            rf=SR(expr).subs(self.x==pr.r_series(self.x,n)+pr.p_series(self.x,n)).subs(t==mu*tau)
            poly=rf.series(mu,n+1).truncate()
            return poly.subs([eq.subs(t==0) for eq in pr.r_solve(n)])
    def p_list(self,n):
        pr=singular_problem(self.x,self.f,self.g)
        return pr.p_series(self.x,n).coefficients(sparse=False)
    def p_odes(self,n):
        pr=singular_problem(self.x,self.f,self.g)
        ode=diff(pr.p_series(self.x,n),tau)-pr.p_series(self.f,n)
        return ode.coefficients(sparse=False)
    def p_ics(self,n):
        pr=singular_problem(self.x,self.f,self.g)
        R=pr.r_list(n)
        P=pr.p_list(n)
        G=[self.g.subs(mu=0)]+[self.g.coefficient(mu^nn) for nn in range(1,n+1)]
        S=pr.r_solve(n)
        ans=[]
        for (r,p,g) in zip(R,P,G):
            ans.append((p==g-r).subs(S).subs([t==0, tau==0]))
        return ans
