%% Experimento Calibración
% Funciona con 2 arduinos: master_tone_ingka.ino y slave_servo.ino
% Testea diferentes fiteos de las asincronías temporales de cada escalon
clear all;
delete(instrfind);

%% Definiciones Experimento

% Numero de bips por trial.
N_stim = 35;
% Repeticiones por condicion.
n=10;
n_calibracion=10;
n_entrenamiento=1;
% Parametros Calibración
out_limit=150;
% Perturbaciones mecanicas
mech_sizes = [33 95]; % tamaï¿½o de la perturbaciï¿½n
servo_ini_pos=57;
% mech_bip_range = [10 13];    % rango de bip
mech_bip_range = [15 18];    % rango de bip

cond_mech = max(size(mech_sizes));

%% 1. Registra datos del sujeto
dir=pwd;
exp='simetrico'; % sufijo del .mat del experimento
filename=['sujetos_' exp '.mat'];
[ sujeto,sujeto_number ] = f_DatosSujeto(dir, exp);





%% Entrenamiento
tent = fopen('temporal_simetrico_ent','w'); % Abre archivos temporales
temp_sizes=[-20 20];
step=1;
conditions='full';
[trial,bad,mech_sizes,temp_sizes] = Loop_central_SIMETRICO(conditions,exp,temp_sizes,step,N_stim,n,n_entrenamiento,n_calibracion,servo_ini_pos,cond_mech,mech_sizes,mech_bip_range);
% Guarda los datos
sujeto(sujeto_number).ent=trial;
sujeto(sujeto_number).ent_mech_sizes=mech_sizes;
sujeto(sujeto_number).ent_temp_sizes=temp_sizes;
clear trial

if  isempty(bad) == 0
    sujeto(sujeto_number).ent_bad=bad;
    clear bad
else
end

disp('Fin del entrenamiento. Presione una tecla para comenzar con el experimento');
disp(' ');
pause()

%% Guarda todos los datos
save(filename,'sujeto')

%% Calibración
tcal = fopen('temporal_simetrico_cal','w');
step=2;
temp_sizes=0;
mech_sizes=mech_sizes(mech_sizes~=57);
conditions='mech';
[trial,bad,mech_sizes,temp_sizes] = Loop_central_SIMETRICO(conditions,exp,temp_sizes,step,N_stim,n,n_entrenamiento,n_calibracion,servo_ini_pos,cond_mech,mech_sizes,mech_bip_range);
% Guarda los datos
sujeto(sujeto_number).cal=trial;
sujeto(sujeto_number).cal_mech_sizes=mech_sizes;
sujeto(sujeto_number).cal_temp_sizes=temp_sizes;
clear trial

if  isempty(bad) == 0
    sujeto(sujeto_number).cal_bad=bad;
    clear bad
else
end

disp('Fin de la primera parte.! Presione una tecla para continuar');
disp(' ');
pause()

%% Experimento
texp = fopen('temporal_simetrico_exp','w');
step=3;
[temp_sizes,esc_sizes] = Parametros_calibracion_exp(sujeto_number,sujeto,out_limit,mech_sizes);
altura=esc_sizes;
[ angulo ] = alt2ang_servo( altura );
mech_sizes=angulo;
conditions='full-iso';
[trial,bad,mech_sizes,temp_sizes] = Loop_central_SIMETRICO(conditions,exp,temp_sizes,step,N_stim,n,n_entrenamiento,n_calibracion,servo_ini_pos,cond_mech,mech_sizes,mech_bip_range);
% Guarda los datos
sujeto(sujeto_number).exp=trial;
sujeto(sujeto_number).exp_mech_sizes=mech_sizes;
sujeto(sujeto_number).exp_temp_sizes=temp_sizes;
clear trial

if  isempty(bad) == 0
    sujeto(sujeto_number).exp_bad=bad;
    clear bad
else
end

disp('Fin del experimento.¡Muchas gracias por participar!');
disp(' ');

%% Guarda los datos en el mat en común
save(filename,'sujeto')
%% Guarda los datos del sujeto en su propia carpeta
mkdir(fullfile('data',num2str(sujeto_number)))
movefile('temporal_simetrico_ent.mat',fullfile('data',num2str(sujeto_number),[num2str(sujeto_number) '_data_simetrico_ent.mat']));
movefile('temporal_simetrico_cal.mat',fullfile('data',num2str(sujeto_number),[num2str(sujeto_number) '_data_simetrico_cal.mat']));
movefile('temporal_simetrico_exp.mat',fullfile('data',num2str(sujeto_number),[num2str(sujeto_number) '_data_simetrico_exp.mat']));