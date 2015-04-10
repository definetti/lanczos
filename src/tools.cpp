#include "tools.h"
#include <stdlib.h>

int Tools::factorial(int n)
{
    int j=n-1;
    if(n==0||n==1)
    {
        return 1;
    }
    do
    {
        n=n*j;
        j--;
    }while (j>0);
    return n;
}

int Tools::binomial(int n, int k)
{
    return factorial(n)/(factorial(k)*factorial(n-k));
}

int Tools::filling(ivec v)
{
    int i,a=0;
    for(i=0;i<(int)v.n_elem;i++)
    {
        if(v(i)!=0)
        {
            a++;
        }
    }
    return a;
}

int Tools::count_double(ivec v, ivec w)
{
    int j,n,q;
    n=0;
    for(j=0;j<n_sites;j++)
    {
        q=v(j)+w(j);
        if(q==2)
        {
            n++;
        }
    }
    return n;
};

int Tools::are_equal(ivec v, ivec w)
{
    int i;
    for(i=0;i<n_sites;i++)
    {
        if(v(i)!=w(i))
        {
            return 0;
        }
    }
    return 1;
}

int Tools::are_neighbours(ivec v, ivec w)
{
    int i,j;

    if (filling(v)!=filling(w)) //initial condition: if they have different electrons throw away
    {
        return 0;
    }

    if(are_equal(v,w)==1) //other condition: if they are equal throw away
    {
        return 0;
    }

    if(v(0)-w(0)!=0&&v(n_sites-1)-w(n_sites-1)==w(0)-v(0)) //PBC hopping case
    {
        for(j=1;j<n_sites-1;j++)
        {
            if(v(j)-w(j)!=0)
            {
               return 0;
            }
        }
        return 1;
    }

    else
    {
        i=0;
        while (v(i)-w(i)==0)
        {
            i++;
        }
        if(v(i)-w(i)!=w(i+1)-v(i+1))
        {
            return 0;
        }
        for(i=i+2;i<n_sites;i++)
        {
            if(v(i)-w(i)!=0)
            {
                return 0;
            }
        }
        return 1;
    }
};

ivec Tools::dec_to_bin(int a, ivec v)
{
    int i,j,q;
    v=zeros<ivec>(n_sites);
    a=a+1;
    i=0;
    q=0;
    for(j=0;j<n_sites;j++)
    {
        if(q<n_electrons&&a>=0)
        {
            if(a<=binomial(n_sites-i-1,n_electrons-q))
            {
                v(i)=0;
            }
            else
            {
                v(i)=1;
                a=a-binomial(n_sites-i-1,n_electrons-q);
                q++;
            }
            i++;
        }
    }
    return v;
}

couple Tools::get_addresses(int i)
{
    int dim = binomial(n_sites,n_electrons);
    couple c;
    c.up= (int)(i / dim);
    c.down=(i % dim);
    return c;
}

int Tools::retrieve_addresses(couple v)
{
    int n;
    int dim=binomial(n_sites,n_electrons);
    n=v.up*dim+v.down;
    return n;
}

void Tools::generate_H() //it is onlt Hup
{
    int i,j,n=0;
    int dim=binomial(n_sites,n_electrons);
    H=zeros<imat>(1,2);
    vector1=zeros<ivec>(n_sites);
    vector2=zeros<ivec>(n_sites);

    n=0;
    for(i=0;i<dim;i++)
    {
        vector1=dec_to_bin(i,vector1);
        for(j=i;j<dim;j++)
        {
            vector2=dec_to_bin(j,vector2);
            if(are_neighbours(vector1,vector2)==1)
            {
                n++;
                H.insert_rows(n,1);
                H(n,0)=i;
                H(n,1)=j;
            }
        }
    }
    H.shed_row(0);
}

void Tools::generate_next()
{
    int i,j,q;
    int dimension=binomial(n_sites,n_electrons)*binomial(n_sites,n_electrons);
    next=zeros<vec>(dimension);

    for(i=0;i<dimension;i++)
    {
        if(curr(i)!=0)
        {
            spins=get_addresses(i);

            next(i)=next(i)+curr(i)*U*count_double(dec_to_bin(spins.up,vector1),dec_to_bin(spins.down,vector2));

            for(j=0;j<(int)H.n_rows;j++) //checks neighbours for the up
            {
                if(H(j,0)==spins.up)
                {
                    spins.up=H(j,1);
                    q=retrieve_addresses(spins);
                    next(q)=next(q)+curr(i)*(-t);
                    spins=get_addresses(i);
                }
                else if(H(j,1)==spins.up)
                {
                    spins.up=H(j,0);
                    q=retrieve_addresses(spins);
                    next(q)=next(q)+curr(i)*(-t);
                    spins=get_addresses(i);
                }
                if(H(j,0)==spins.down) //checks the neighbour for the down
                {
                    spins.down=H(j,1);
                    q=retrieve_addresses(spins);
                    next(q)=next(q)+curr(i)*(-t);
                    spins=get_addresses(i);
                }
                else if(H(j,1)==spins.down)
                {
                    spins.down=H(j,0);
                    q=retrieve_addresses(spins);
                    next(q)=next(q)+curr(i)*(-t);
                    spins=get_addresses(i);
                }
            }
        }
    }
}

vec Tools::build_lancmat()
{
    int counter;
    double alpha=0,beta=0;
    int dimension=binomial(n_sites,n_electrons)*binomial(n_sites,n_electrons);
    lancmat=zeros<mat>(n_iter,n_iter);
    vector1=zeros<ivec>(n_sites);
    vector2=zeros<ivec>(n_sites);
    prev=zeros<vec>(dimension);
    curr=zeros<vec>(dimension);
    next=zeros<vec>(dimension);
    vec eigval(n_iter);

    curr(0)=1;

    for(counter=0;counter<n_iter-1;counter++)
    {

        generate_next();
        alpha=dot(next,curr);
        next=next-alpha*curr-beta*prev;
        beta=norm(next);
        lancmat(counter,counter)=alpha;
        lancmat(counter+1,counter)=beta;
        lancmat(counter,counter+1)=beta;
        prev=curr;
        curr=next/beta;
    }

    generate_next();
    lancmat(n_iter-1,n_iter-1)=dot(next,curr);
    eigval=eig_sym(lancmat);

    return eigval;
};
