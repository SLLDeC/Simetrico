clear all
close all
load('sujetos_simetrico.mat')

%% Par�metros----------------------------------------------
pre_baseline_bips = 5;        % N° de bips usados para calcular el prebaseline
pos_baseline_bips = 5;        % N° de bips usados para calcular el posbaseline
out_limit = 150;              % Asincronia maxima permitida
%%----------------------------------------------------------
n_rep=10; % num. de repeticiones por condicion
c=length(sujeto); % num. de sujetos
max_trial=[-17:23]; % coordenadas del trial
ASYN=nan(n_rep,length(max_trial),c);
ASYN_BS=nan(n_rep,length(max_trial),c);

for s=1:length(sujeto)
    
    ASYN=nan(n_rep,length(max_trial),c);
    ASYN_BS=nan(n_rep,length(max_trial),c);
    
    N=length(sujeto(s).exp);
    
    cposITI=1;cpreITI=1;
    
    % Arma los vectores para cada trial
    for m=1:N
        mech_sizes(m)=sujeto(s).exp(m).mech_size; % vector de pert. mech
        pert_bip(m)=sujeto(s).exp(m).mech_bip;
        pert_indx(m)=find(sujeto(s).exp(m).asyn(2,:)==pert_bip(m));
        min_bip(m,sujeto(s).exp(m).mech_size)=min(sujeto(s).exp(m).asyn(2,:)-pert_bip(m));
        max_bip(m,mech_sizes(m))=max(sujeto(s).exp(m).asyn(2,:)-pert_bip(m));
        pre_baseline(1,m)=mean(sujeto(s).exp(m).asyn(1,pert_indx(m)-pre_baseline_bips:pert_indx(m)));
        pos_baseline(m)=mean(sujeto(s).exp(m).asyn(1,end-pos_baseline_bips:end));
        
        if isempty(find(abs(sujeto(s).exp(m).asyn(1,:)-pre_baseline(m))>=out_limit))==1
            out(m)=0;
        else
            out(m)=1;
        end
    end
    
    steps(1,:)=sort(repmat(sort(sujeto(s).exp_mech_sizes),1,3));
    steps(2,:)=repmat(sort(sujeto(s).exp_temp_sizes),1,3);
    
    Ctemp=steps(2,:)==0;
    Cmech=steps(1,:)==57;
    Cond=Cmech & Ctemp;
    steps(:,Cond)=[];
    
    contador=ones(1,length(steps));
    
    
    for n=1:length(steps)
        min_bip_cond(n)=max(nonzeros(min_bip(:,steps(1,n))));
        max_bip_cond(n)=min(nonzeros(max_bip(:,steps(1,n))));
    end
    
    for t=1:N
        
        sujeto(s).exp(t).out=out(t);
        
        if out(t)==0 % Si el trial no tiene outlier se calcula el pre y posITI
            for i=1:length(sujeto(s).exp(t).resp)-1
                resp_n=sujeto(s).exp(t).resp(1,i);
                resp_n1=sujeto(s).exp(t).resp(1,i+1);
                if i<pert_bip(t)
                    preITI(cpreITI)=resp_n1-resp_n;
                    cpreITI=cpreITI+1;
                elseif i>=pert_bip(t)+5 %el posITI se calcula luego del "overshoot"
                    posITI(cposITI)=resp_n1-resp_n;
                    cposITI=cposITI+1;
                end
            end
            meanpreITI(t)=mean(preITI);
            meanposITI(t)=mean(posITI);
        else % Si el trial tiene un outlier el ITI se computa como NaN
            meanpreITI(t)=nan;
            meanposITI(t)=nan;
        end
        
        sizes=[sujeto(s).exp(t).mech_size sujeto(s).exp(t).temp_size];
        for h=1:length(steps)
            if isequal(sizes,steps(:,h)')==1
                c=h;
                contador(h)=contador(h)+1;
                f=contador(h);
            end
        end
        
        idx=sujeto(s).exp(t).asyn(2,1)-pert_bip(t)+18;
        sujeto(s).condicion(c).mech_size=steps(c);
        if out(t)==0
            ASYN(f,idx:length(sujeto(s).exp(t).asyn(1,:))+idx-1,c)=sujeto(s).exp(t).asyn(1,:);
            ASYN_BS(f,idx:length(sujeto(s).exp(t).asyn(1,:))+idx-1,c)=sujeto(s).exp(t).asyn(1,:)-pre_baseline(t);
            PRE_BS(f,c)=pre_baseline(t);
            POS_BS(f,c)=pos_baseline(t);
            PRE_ITI(f,c)=meanpreITI(t);
            POS_ITI(f,c)=meanposITI(t);
        else
            ASYN(f,idx:length(sujeto(s).exp(t).asyn(1,:))+idx-1,c)=nan(1,length([idx:length(sujeto(s).exp(t).asyn(1,:))+idx-1]));
            ASYN_BS(f,idx:length(sujeto(s).exp(t).asyn(1,:))+idx-1,c)=nan(1,length([idx:length(sujeto(s).exp(t).asyn(1,:))+idx-1]));
            PRE_BS(f,c)=nan;
            POS_BS(f,c)=nan;
            PRE_ITI(f,c)=nan;
            POS_ITI(f,c)=nan;
        end
        
    end
    
    
    for c=1:length(sujeto(s).condicion)
        
        sujeto(s).condicion(c).serie_prom(1,:)=nanmean(ASYN(:,:,c));
        sujeto(s).condicion(c).serie_prom(2,:)=max_trial;
        sujeto(s).condicion(c).serie_prom_std(1,:)=nanstd(ASYN(:,:,c));
        sujeto(s).condicion(c).serie_prom_std(2,:)=max_trial;
        sujeto(s).condicion(c).serie_prom_bs(1,:)=nanmean(ASYN_BS(:,:,c));
        sujeto(s).condicion(c).serie_prom_bs(2,:)=max_trial;
        sujeto(s).condicion(c).serie_prom_bs_std(1,:)=nanstd(ASYN_BS(:,:,c));
        sujeto(s).condicion(c).serie_prom_bs_std(2,:)=max_trial;
        sujeto(s).condicion(c).mean_prebaseline=nanmean(PRE_BS(:,c));
        sujeto(s).condicion(c).std_prebaseline=nanstd(PRE_BS(:,c));
        sujeto(s).condicion(c).mean_posbaseline=nanmean(POS_BS(:,c));
        sujeto(s).condicion(c).std_posbaseline=nanstd(POS_BS(:,c));
        sujeto(s).condicion(c).mean_preITI=nanmean(PRE_ITI(:,c));
        sujeto(s).condicion(c).std_preITI=nanstd(PRE_ITI(:,c));
        sujeto(s).condicion(c).mean_posITI=nanmean(POS_ITI(:,c));
        sujeto(s).condicion(c).std_posITI=nanstd(POS_ITI(:,c));
        
        %     xlswrite('ASYN_2.xls',ASYN(:,:,2))
        
        
        %     sujeto(s).condicion(1).n=c1;
        
        
    end
    
    
    save('sujetos_alturas.mat','sujeto')
    
    %     %% Series temporales de cada sujeto
    colores=['g' 'r' 'k' 'm' 'c' 'b' 'y'];
    figure(s)
    for p=1:length(sujeto(s).condicion)
        
        [lineprops] = fg_Colores_presentacion(p);
        lineprops=cell2mat(lineprops.col(1));
        title(['sujeto ' num2str(s)])
        %         color=colores(p);
        %         plot(sujeto(s).condicion(p).serie_prom(2,:),sujeto(s).condicion(p).serie_prom(1,:)-sujeto(s).condicion(p).mean_prebaseline,'.-','Color',color)
        hold on
        plot(sujeto(s).condicion(p).serie_prom(2,:),sujeto(s).condicion(p).serie_prom_bs(1,:),'Color',lineprops)
        hold all
    end
    legend('-M-T','-M','-M+T','-T','+T','+M-T','+M','+M+T')
    %
    %
end

%% Efecto de la perturbacion espacial en la asincronia
figure()


for s=1:length(sujeto)
    
    x_mech=sort(repmat(sujeto(s).exp_mech_sizes,1,3));
    x_mech(4:6)=[];
    
    color=colores(s);
    Legend{s}=strcat('sujeto', num2str(s));
    for p=1:length(sujeto(s).condicion)
        title('pert. espacial')
        ylabel('asincron�a [ms]')
        xlabel('altura escal�n [cm]')
        data(p)=sujeto(s).condicion(p).serie_prom_bs(1,(sujeto(s).condicion(p).serie_prom(2,:)==1));
        error(p)=sujeto(s).condicion(p).serie_prom_bs_std(1,(sujeto(s).condicion(p).serie_prom(2,:)==1));
    end
    
    data(4:5)=[];
    error(4:5)=[];
    errorbar(x_mech,data,error,'.-','LineWidth',2,'MarkerSize',5,'Color',color)
    hold on
    %     xlim([-2 2])
     title('Asyn')
end
legend(Legend)

%% Efecto de la perturbacion espacial en deltaNMA
figure()
for s=1:length(sujeto)
    [lineprops] = fg_Colores_presentacion(s);
    lineprops=cell2mat(lineprops.col(1));
    Legend{s}=strcat('sujeto', num2str(s));
    color=colores(s);
    for p=1:length(sujeto(s).condicion)
        title('deltaNMA')
        ylabel('posNMA-preNMA [ms]')
        xlabel('altura escalón [cm]')
        data(p)=sujeto(s).condicion(p).mean_prebaseline-sujeto(s).condicion(p).mean_posbaseline;
    end
    data(4:5)=[];
    error(4:5)=[];
    plot(x_mech,data,'.-','LineWidth',2,'MarkerSize',5,'Color',color)
    hold on
    %     xlim([-2 2])
     title('deltaNMA')
end
legend(Legend)

%% Efecto de la perturbacion espacial en deltaITI
figure()
for s=1:length(sujeto)
    color=colores(s);
    Legend{s}=strcat('sujeto', num2str(s));
    for p=1:length(sujeto(s).condicion)
        
        ylabel('posITI-preITI [ms]')
        xlabel('altura escal�n [cm]')
        data(p)=sujeto(s).condicion(p).mean_preITI-sujeto(s).condicion(p).mean_posITI;
    end
    data(4:5)=[];
    plot(x_mech,data,'.-','LineWidth',2,'MarkerSize',5,'Color',color)
    
    hold on
    %     xlim([-2 2])
    title('deltaITI')
end
legend(Legend)

save('sujetos_simetrico.mat','sujeto')