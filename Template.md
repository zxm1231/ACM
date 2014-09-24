##连通性
###SCC	
	void tarjan(int u)
	{
	    dfn[u] = low[u] = dindex++;
	    instack[u] = 1;
	    st[++tail] = u;
	    for(int j = head[u] ; ~j ; j = e[j].next){
	        int v = e[j].v;
	        if (!dfn[v])
	        {
	            tarjan(v);
	            if (low[u] > low[v])low[u] = low[v];
	        }
	        else if (instack[v] && dfn[v] < low[u])
	            low[u] = dfn[v];
	    }
	    if (dfn[u] == low[u]){
	        int i;
	        do{
	            i = st[tail--];
	            instack[i] = 0;
	            cmp[i] = bcnt;
	        }
	        while (i!=u);
	        bcnt++;
	    }
	}
	void scc()
	{
	    bcnt = dindex = 1;
	    tail = 0;
	    clr(dfn , 0);
	    clr(low , 0);
	    clr(instack , 0);
	    clr(st , 0);
	    clr(cmp , 0);
	    for (int i = 1 ; i <= N ; i++)
	        if (!dfn[i])tarjan(i);
	}

