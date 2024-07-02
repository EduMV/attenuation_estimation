function ps_data = bf_planewave(sensor_data_kwave, fs)

data = sensor_data_kwave';

dx = 0.3e-3;
x=(1:128)*dx; %%% posiciones en x de la imagen
x = x-mean(x);
c = 1540;
%fs = 1/1.168831168831169e-08;
tx_positions = x;


for j=1:length(x) %% Iterar en x
    
    %h1 = waitbar(j/length(x),h1,['Generando líneas... ' num2str(j) '/' num2str(length(x))]);
    
    for i=1:size(data,1) %% Iterar en z
        
        %%% Posición del punto en x,z
        
        x_point=x(j);
        z_point=c*(i/fs)/2;
        z_rx=0;
        
        
        %%% Calcular retardos
        ang = 0;
        d_ec=z_point*cos(ang)+x_point*sin(ang);
        tau=(d_ec+sqrt(z_point^2+(x_point-tx_positions).^2))/c;
        
        %%% Calcular muestras (no se usa interpolación)
        
        muestra=round((tau)*fs);
        muestra(muestra<=0)=1;
        muestra(muestra>size(data,1))=1;
        ind=muestra+size(data,1)*(0:size(muestra,2)-1); %Indexación linear
        values=data(ind);
        
        %%% Sumar valores solo de la apertura
        F = 3;
        a=z_point/(2*F);
        liminf=x_point-a;
        limsup=x_point+a;
        apodization=(tx_positions<limsup & tx_positions>liminf);
        
        ps_data(i,j,1)=sum(apodization*values');
        
    end
    %keyboard

end