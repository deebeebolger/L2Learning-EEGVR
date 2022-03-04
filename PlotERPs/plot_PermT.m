function[t, pvals_corr,X] = plot_PermT(Data,time,lim,green_on)
% *********************************************************************
% Date: Octobre 2014            Programm� par: D. Bolger
% Function to carry out a permutation test with correction fdr entre
% deux vecteur de donn�es (le grande moyen pour deux conditions 
% sujets).
% La fonction cherche les pointes temporelles o� les p-values sont <=0.05
% et les pointes temporelles o� les p-values sont <0.05 et <0.06. Ces
% pointes temporelles sont marqu�es en rouge et verts, respectivement. 
% Data:  matrice (Time X 2) 
% time: time vector
% lim: d�termine � quelle niveau dans le plot les p-values sont pr�sent�.
% On peut le d�finir pour que la repr�sentation des p-values ne g�ne pas la
% reste du plot, e.g. 4.
% Output:
% pvals_corr = vecteur de permutation p-values (m�me taille que le vecteur time). 
%**************************************************************************************

dt = diff(time);
[t,df,pvals]=statcond(Data,'mode','perm','naccu',5000);  %calculated with a default number of 1000 random partitions
[pvals_corr, pvals_masked] = fdr(pvals);   % fdr correction of permutation p-values
assignin('base','pvals_corr',pvals_corr);
assignin('base','pvals',pvals);
%pvals_corr(ti)=1;
xi=(pvals_corr<=0.05);
if ~isempty(xi)
    xi_diff = diff(xi);
    xi1 = find(xi_diff == 1);
    xi2 = find(xi_diff == -1);
    if length(xi1) > length(xi2)
        xi2 = cat(1,xi2,length(xi));
    end
    X=zeros(length(time),1);
    for cntr = 1:length(xi1)
        d = (xi2(cntr) - xi1(cntr))*dt(1);
        if d<10
            X(xi1:xi2) = 0;
        elseif d>=10
            X(xi1(cntr):xi2(cntr)) = lim;
        end
        
    end
end

xi2=find(pvals_corr>0.05 & pvals_corr<=0.06);

X2=zeros(length(time),1);

i=find(X==0);
X(i)=NaN;

X2(xi2)=lim;
i1=find(X2 ==0);
X2(i1)=NaN;

p5=plot(gca,time,X,'ro');

if strcmp(green_on,'yes')
    hold on
    p6=plot(gca,time,X2,'go');
end

end

