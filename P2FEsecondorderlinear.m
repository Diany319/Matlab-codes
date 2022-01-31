function x=P2FEsecondorderlinear(Mg,u0,u1,C)
    %Ce pragramme calcule la solution approchee en 1D en elements finis P2
    % du probleme de Dirichlet non homogene -Cu"+u =f sur (0,1) 
    %avec les conditions aux limites u(0)= u0, u(1)=u1
    % Mg est une matrice 1 ligne (N+1) colonnes constituee des noeuds geometriques dans ]0,1[.
    % C est une constante
     
    mG=[-sqrt(3/5), sqrt(3/5), 0 ];% la matrice des points d'integration
    w=[5/9, 5/9, 8/9];  % la matrice  des poids d'integration

    [Phi,dPhi]= lagrangeinterpo(mG);
    
    %Nombre d'elements du maillage
    Nel=length(Mg)-1;
    
    %Vecteur des longueurs des elements
    h=Mg(2:Nel+1)-Mg(1:Nel);
    
    %Tableau de connectivite
    Connec=zeros(Nel,3);
    for K=1:Nel
	  Connec(K,1:3)= [2*K-1 2*K+1 2*K];
    end
    
    %Tableau de numerotation
    Numer=zeros(2*Nel+1,1);
    Numer(1)=2*Nel;
    Numer(2*Nel+1)=2*Nel+1;
    
    Numer(2:2*Nel)=1:2*Nel-1;
 
   %Tableau d'adressage
    adres=zeros(Nel, 3);
    for K=1:Nel  
      for i=1:3
         adres(K,i)=Numer(Connec(K,i));
      end
    end

    % assemblage des matrices locales 
     A= zeros(2*Nel+1,2*Nel+1); F= zeros(2*Nel+1,1);
     
    for K=1:Nel % boucle sur les elements
       for i=1:3 % boucle sur les lignes
	    
            floc(i)= (h(K)/2)*sum(w(1:3).*(f(Mg(K+1)+Mg(K)+ h(K)*mG(1:3)).*Phi(i,1:3)))-1*(h(1)/2)*sum(w(1:3).*(u0*Phi(1,1:3).*Phi(i,1:3)+u1*Phi(2,1:3).*Phi(i,1:3)))-((2*C)/h(Nel))*sum(w(1:3).*(u0*dPhi(1,1:3).*dPhi(i,1:3)+u1*dPhi(2,1:3).*dPhi(i,1:3)));
        
            for j=1:3 % boucle sur les colonnes
                Aloc(i,j)=h(K)/2*sum(w(1:3).*(Phi(i,1:3).*Phi(j,1:3)))+((2*C)/h(K))*sum(w(1:3).*(dPhi(i,1:3).*dPhi(j,1:3)));
            end
      end
         
	 %assemblage dans la matrice globale et du second menbre
         for i=1:3 
            for j=1:3
                A(adres(K,i),adres(K,j)) = A(adres(K,i),adres(K,j)) + Aloc(i,j);
            end
         F(adres(K,i))=F(adres(K,i))+ floc(i);
         end
         
    end 
  
     %solution du syst�me global

     S = A(1:2*Nel-1,1:2*Nel-1)\F(1:2*Nel-1)
     u(2:2*Nel) = S;
     u(1)= u0;
     u(2*Nel+1)=u1;

    %Graphe de la solution passant par les noeuds de calcul

    % hold on;
     xx = zeros(2*Nel+1,1); 
     for K =1:Nel
       xx(2*Nel+1) = Mg(Nel+1);
       xx(2*K-1) = Mg(K);
       xx(2*K)   = (Mg(K+1)+Mg(K))/2;
     end
     
     %f=((24-35*exp(1))*exp(xx)+(-24*exp(2)+35*exp(1))*exp(-xx))/(1-exp(2))+xx.^4-12*xx.^2-24
     %plot(xx,u, 'b',xx,f, 'g')
     plot(xx,u, 'b')
     hold on
    plot(xx,uExact(xx),'r')
     hold on;
     legend('Sol num', 'Sol exact');
     %axis([0,5,-15,5]);
     ff=getframe;
    hold off;
 end   
 
% definition des sous-fonctions 

    function [Psi,dPsi]= lagrangeinterpo(mG)
    %Psi et dPsi sont des matrices des valeurs des polynomes 
    % de Lagrance sur l'element de reference [-1,1] 
     % aux points  de la quadrature de Gauss.
       nndG=length(mG);
       Psi=zeros(3,nndG); dPsi =zeros(3,nndG);
       for iG=1:nndG
           xi=mG(iG);
	       Psi(1,iG)=xi*(xi-1)/2;
	       Psi(2,iG)=xi*(xi+1)/2;
	       Psi(3,iG)=1-xi^2;
	       dPsi(1,iG)= (2*xi-1)/2;
           dPsi(2,iG)= (2*xi+1)/2;
           dPsi(3,iG)= -2*xi;
       end
    end 


   function f=f(x)
     %f=x.^4;
    % f = x;
    f=0.0
   end
