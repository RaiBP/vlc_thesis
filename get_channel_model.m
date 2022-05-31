function [channel_impulse_response, fs] = get_channel_model(tx_position, rx_position)
    %posición Ti
    x_i=tx_position(1);
    y_i=tx_position(2);
    z_i=tx_position(3);

    %posición Rj1    
    x_j=rx_position(1);
    y_j=rx_position(2);
    z_j=rx_position(3);

    %posición del área reflectiva (pared) usar para puntos especificos de
    %elementos reflectantes

    % % first reflection_point 
    x_w1=2;
    y_w1=2;
    z_w1=3;
    % % second reflection_point
    x_w2=1.5;
    y_w2=2.3;
    z_w2=4.0;
    % % third reflection_point
    x_w=5;
    y_w=2.9;
    z_w=3.5;

    % parámetros de la simulación %

    Aw=1  ;                                     %área del elemento reflectante
    pw=0.6 ;                                    %coeficiente de reflexión del área reflectiva
    ang_rad = 60;                               %semi-ángulo de mitad de potencia del LED
    m= -log(2)/log(abs(cos(ang_rad*pi/180)));   %número Lambertiano
    Ap=0.0001;                                  %área física del receptor (1 cm^2)
    eta=1.5;                                    %indice de refracción del PD
    fov = 70;                                    %field of view
    
    %parámetros de los ángulos de inclinacin y rotacin%
    beta_i=45;                                  %ángulo de inclinación del LED con respecto al eje z
    alpha_i=45;                                 %ángulo de rotación del LED con respecto al eje x
    
    beta_j=45;                                  %ángulo de inclinación del PD con respecto al eje z
    alpha_j=45;                                 %ángulo de rotación del PD con respecto al eje x
    
    % wall roation angles  para angulo especifico                   
    alpha_w=20;                               
    beta_w=70;
    
    gv=[20 60];                                 %uniform distribution probability  
    fv=[0 5];                                   %uniform distribution probability   
    W=6;                                        % width of obstacle
    H=5;                                        % hight of obstacle
    X=6;                                        %Lenght of mine
    Y=3;                                        %height of mine
    es=5;                                       % epsilon 
    t=5*10^-9;                                  % value of time
    gymma=0.017;                                %Reflection coefficient
    g=0.72;                                     %Responsivity
    f=0.5;                                      %Refractive index
    kr=[0.1 0.01];                              %uniform distribution probability
    km=[0 10];                                  %uniform distribution probability
    ks=kr+km;                                   %uniform distribution probability
    N=70;                                       %Number de scatters

    p=0.1;                                      % parameter is used in calculation of Gn below equation 3.21 
    
    c=3*10^8 ;                                  % speed of light

    sampling_time=0.25e-9; %(0.25 nano segundo de muestreo, se disminuye para mayor numero de muestras)
    fs = 1/sampling_time;
    time=0:sampling_time:35e-9; % vector de tiempo donde se observará el CIR
    time=round(time,12);    %redondeo del tiempo a 12 cifras significativas 
    t_rise=0.5e-9; %time rise subida del PD
    t_fall=1e-9;    %time bajada del PD
    h_led=10*(exp(-time/t_fall)-exp(-time/(t_rise))); %respuesta al impulso del PD

    %------Angulos de incidencia e irradiancia entre los enlaces------%

    %angulo de incidencia del enlace entre LED y PD para 5 posiciones del PD

    %1st position of receiver
    incidencia_radian1=incline(x_i,y_i,z_i,x_j,y_j,z_j,alpha_j,beta_j);
    incidencia1=rad2deg(incidencia_radian1) ; %conversion de radianes a grados


    %angulo de irradiancia del enlace entre LED y PD para 5 posiciones del PD
    %1st position of receiver
    irradiancia_radian1=rotacion(x_j,y_j,z_j,x_i,y_i,z_i,alpha_i,beta_i);
    irradiancia1=rad2deg(irradiancia_radian1); %conversion de radianes a grados

    %angulo de irradiancia del enlace entre LED y pared 
    irradiancialw_radian=rotacion(x_w,y_w,z_w,x_i,y_i,z_i,alpha_i,beta_i);
    irradiancialw=rad2deg(irradiancialw_radian) ; %conversion de radianes a grados

    %angulo de incidencia del enlace entre LED y pared
    incidencialw_radian=incline(x_i,y_i,z_i,x_w,y_w,z_w,alpha_w,beta_w);
    incidencialw=rad2deg(incidencialw_radian);  %conversion de radianes a grados

    %angulo de irradiancia del enlace entre pared y PD
    irradianciaw_radian=rotacion(x_j,y_j,z_j,x_w,y_w,z_w,alpha_w,beta_w);
    irradianciaw=rad2deg(irradianciaw_radian);  %conversion from radian to degree

    %angle de incidencia del enlace de pared a PD
    %1st position of receiver
    incidenciaw_radian11=incline(x_w1,y_w1,z_w1,x_j,y_j,z_j,alpha_j,beta_j);
    incidenciaw11=rad2deg(incidenciaw_radian11);%conversion de radianes a grados

    incidenciaw_radian21=incline(x_w2,y_w2,z_w2,x_j,y_j,z_j,alpha_j,beta_j);
    incidenciaw21=rad2deg(incidenciaw_radian21);

    incidenciaw_radian31=incline(x_w,y_w,z_w,x_j,y_j,z_j,alpha_j,beta_j);
    incidenciaw31=rad2deg(incidenciaw_radian31);

    % Calculo de la distancia euclidiana entre LED y PD para 5 puntos 
    D1 = DIST(x_i,y_i,z_i,x_j,y_j,z_j);
    
    %-------- Calculo del CIR para la componente LoS------%
    h_vector = zeros(1,length(time)); %vector donde se guarda la CIR del LoS

    %Funci�n que calcula el CIR y el delay del LoS
    [m_HLoS1,dl1] = HLoS_direct(x_i,y_i,z_i,x_j,y_j,z_j,Ap,eta,alpha_i,alpha_j,beta_i,beta_j,incidencia_radian1,incidencia1,m,fov,gv,fv,W,H,X,Y,t,es,c);
    ind = round(dl1,9);      %redondeo del delay a 9 cifras significativas
    index = find(ind==time); %encuentra el indice cuando el delay es igual a la función de tiempo
    h_vector(index) = h_vector(index) + m_HLoS1;  %coloca el CIR en el indice de delay correspondiente
    %Convolucion entre CIR del LoS y el CIR del PD
    Final_responsel = conv(h_led, h_vector);
    
    
    %-------- Calculo del CIR para la componente nLoS de una pared total fija en eje y------%

    a=1; %indice para ir guardando componentes nLoS, delay, angulo de incidencia, angulo de inclinación y angulo de rotación, 
    lx=6; ly=3; lz=5;  %Area del escenario Minero
    Nx=lx*3; Ny=ly*3; Nz=round(lz*3); %numero de grid de cada superficie
    dA=lz*ly/(Ny*Nz); % area de cada grid
    x=0:lx/Nx:lx; %vector que cubre todo el grid en x
    y=0:ly/Ny:ly;  %vector que cubre todo el grid en y
    z=0:lz/Nz:lz;  %vector que cubre todo el grid en z
    HnLos=zeros(1,length(time));  %vector donde se guarda la CIR del nLoS
    tem=zeros(1,length(time));    %vector donde se guarda el delay del nLoS
    incide=zeros(1,length(time)); %vector donde se el angulo de indicencia CIR del nLoS
    for kk=1:Nx+1  
        for ll=1:Nz+1
            WP1 = [x(kk) 0.2 z(ll)];    % ubicación de la pared fija en y se mueven en x y z
            r(a) = randi([0 90],1,1); % angulo de rotacion del elemento reflectante aleatorio uniforme
            s(a) = randi([0 90],1,1); % angulo de inclinación del elemento reflectante aleatorio uniforme
            incidenciaw_pru = incline(x(kk),0.2,z(ll),x_j,y_j,z_j,alpha_j,beta_j); %angulo de radiancia del LED al elemento reflectante
            incidenciawpru = rad2deg(incidenciaw_pru);%conversion de radianes a grados
            %Funci�n que calcula el CIR y el delay del nLoS
            [m_HnLoSp,dnp] = HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x(kk),0.2,z(ll),dA/70,pw,alpha_i,alpha_j,r(a),beta_i,beta_j,s(a),Ap,incidenciawpru,incidenciaw_pru,eta,m,fov,gv,fv,W,H,X,Y,t,es,c)  ;
            incide(a) = incidenciawpru; %guarda angulo de incidencia
            HnLos(a) = m_HnLoSp; %guarda CIR del nLoS
            tem(a) = dnp; %guarda delay del nLoS
            a = a + 1;
        end 
    end

    h_vector3 = zeros(1,length(time)); %vector donde se guarda el CIR del nLoS dependiendo del delay
    tem = round(tem,9); %redondeo del vector de delays
    indno = zeros(1,length(time));   %vector donde se guardara el indice del delay del nLoS = tiempo
    indexno = 1;
    Final_responsen = zeros(1, length(conv(h_led, h_vector3))); %vector donde se guarda la conv entre CIR del LED y CIR nLoS
    for i=1:length(time)
        var = tem(i);  %guarda en una variable el delay
        var = round(var,9); %redonde el delay a 9 cifras significativas
        indexno = find(time==var); %encuentra el indice cuando el delay es igual al tiempo
        indno(i) = indexno;   %guarda ese indice en el vector de indices
        h_vector3(indexno) = h_vector3(indexno) + HnLos(i); % coloca el CIR del nlos en el indice correspondiente
    end
    Final_responsen=conv(h_led,h_vector3);

    %-------- Calculo del CIR para la componente nLoS de una pared total fija en eje y------%

    b=1; %indice para ir guardando componentes nLoS, delay, angulo de incidencia, angulo de inclinación y angulo de rotación, 
    HnLos1=zeros(1,length(time));  %vector donde se guarda la CIR del nLoS
    tem1=zeros(1,length(time));    %vector donde se guarda el delay del nLoS
    incide1=zeros(1,length(time)); %vector donde se el angulo de indicencia CIR del nLoS
    for kk=1:Nx+1  
        for ll=1:Nz+1
            WP2 = [x(kk) 3 z(ll)];    % ubicación de la pared fija en y se mueven en x y z
            r1(b) = randi([0 90],1,1); % angulo de rotacion del elemento reflectante aleatorio uniforme
            s1(b) = randi([0 90],1,1); % angulo de inclinación del elemento reflectante aleatorio uniforme
            incidenciaw_pru1 = incline(x(kk),3,z(ll),x_j,y_j,z_j,alpha_j,beta_j); %angulo de radiancia del LED al elemento reflectante
            incidenciawpru1 = rad2deg(incidenciaw_pru1);%conversion de radianes a grados
            %Función que calcula el CIR y el delay del nLoS
            [m_HnLoSp1,dnp1] = HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x(kk),3,z(ll),dA/70,pw,alpha_i,alpha_j,r1(b),beta_i,beta_j,s1(b),Ap,incidenciawpru1,incidenciaw_pru1,eta,m,fov,gv,fv,W,H,X,Y,t,es,c)  ;
            incide1(b) = incidenciawpru1; %guarda angulo de incidencia
            HnLos1(b) = m_HnLoSp1; %guarda CIR del nLoS
            tem1(b) = dnp1; %guarda delay del nLoS
            b = b + 1;
        end 
    end

    h_vector1 = zeros(1,length(time)); %vector donde se guarda el CIR del nLoS dependiendo del delay
    tem1 = round(tem1,9); %redondeo del vector de delays
    indno1 = zeros(1,length(time));   %vector donde se guardara el indice del delay del nLoS = tiempo
    indexno1 = 1;
    Final_responsen1 = zeros(1,length(conv(h_led,h_vector1))); %vector donde se guarda la conv entre CIR del LED y CIR nLoS
    for i=1:length(time)
        var1 = tem1(i);  %guarda en una variable el delay
        var1 = round(var1,9); %redonde el delay a 9 cifras significativas
        indexno1 = find(time==var1); %encuentra el indice cuando el delay es igual al tiempo
        indno1(i) = indexno1;   %guarda ese indice en el vector de indices
        h_vector1(indexno1) = h_vector1(indexno1) + HnLos1(i); % coloca el CIR del nlos en el indice correspondiente
    end
    Final_responsen1 = conv(h_led,h_vector1);

    %-------- Calculo del CIR para la componente escattering-----%
    Hsca = zeros(1,length(time));
    h_vector2 = zeros(1,length(time));
    temsca = zeros(1,length(time));

    %Función que calcula el CIR y el delay del scatter
    [m_Hscat1,ds1] = H_scater(x_i,y_i,z_i,x_j,y_j,z_j,Ap,m,f,g,gymma,kr,km,ks,p,N,irradiancia1,c,alpha_i,beta_i);

    indexs=1;    
    for i=1:length(m_Hscat1)
        Hsca(i) = m_Hscat1(i);
        temsca(i)= ds1(i);   
    end

    for i=1:length(time)
        vars = temsca(i);
        vars = round(vars, 9);
        indexs = find(time==vars);
        h_vector2(indexs) = h_vector2(indexs)+Hsca(i);
    end

    Final_responses = conv(h_led,h_vector2);
    channel_impulse_response = Final_responsel+Final_responsen+Final_responsen1+Final_responses;
end

%----------Funciones usadas------------%
 
% function to plot the graph for multiple impulse
function multi_pulse_plot = graph_draw(dm_total,m_total_Hscat,t1,D1,N)
    for i=1:N
        y1 = dirac(t1-dm_total(i));
        multi_pulse_plot = plottt(y1,t1,m_total_Hscat(i),D1);
    end
end


% function to calculate the incidence angle
function ang_inc = incline(x1,y1,z1,x2,y2,z2,alpha,beta)
    v=[x1-x2,y1-y2,z1-z2];
    Ntilt=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),cosd(beta)];
    d_p=dot_product(v,Ntilt);
    d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    ang_inc=acos(d_p/d);
end


%funci�n para calcular el angulo de irradiancia%
function ang_incidencia = rotacion(x1,y1,z1,x2,y2,z2,alpha,beta)
    v=[x1-x2,y1-y2,z1-z2];
    Ntilt=[cosd(alpha)*sind(beta),sind(alpha)*sind(beta),-cosd(beta)];
    d_p=dot_product(v,Ntilt);
    d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    ang_incidencia=acos(d_p/d);
end


%funcion to calculate the HnLos for three reflection points
function [HnLoS_Total,delta_t] = HnLos_calculation_total(x_i,y_i,z_i,x_j,y_j,z_j,x_w1,y_w1,z_w1,x_w2,y_w2,z_w2,x_w3,y_w3,z_w3,Aw,pw,alpha_i,alpha_j,alpha_w1,alpha_w2,alpha_w3,beta_i,beta_j,beta_w1,beta_w2,beta_w3,Ap,inc1,inc_r1,inc2,inc_r2,inc3,inc_r3,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);  
    [m_HnLoS1,dm1]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w3,y_w3,z_w3,Aw,pw,alpha_i,alpha_j,alpha_w3,beta_i,beta_j,beta_w3,Ap,inc1,inc_r1,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);
    [m_HnLoS2,dm2]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w1,y_w1,z_w1,Aw,pw,alpha_i,alpha_j,alpha_w1,beta_i,beta_j,beta_w1,Ap,inc2,inc_r2,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);
    [m_HnLoS3,dm3]=HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w2,y_w2,z_w2,Aw,pw,alpha_i,alpha_j,alpha_w2,beta_i,beta_j,beta_w2,Ap,inc3,inc_r3,eta,m,fov,gv,fv,W,H,X,Y,t,es,c);
    HnLoS_Total=[m_HnLoS1,m_HnLoS2,m_HnLoS3];
    delta_t=[dm1,dm2,dm3];
end


%funcion que calcula el HNLoS
function [m_HnLoS,dm] = HnLos_calculation(x_i,y_i,z_i,x_j,y_j,z_j,x_w,y_w,z_w,Aw,pw,alpha_i,alpha_j,alpha_w,beta_i,beta_j,beta_w,Ap,inc,inc_r,eta,m,fov,gv,fv,W,H,X,Y,t,es,c)
    dv_iw= dv(x_i,y_i,x_w,y_w,fv);
    sv_iw= sv(x_i,y_i,z_i,x_w,y_w,z_w,fv);
    Piw=P_expt(gv,fv,W,H,X,Y,t,es, dv_iw,sv_iw);
    
    dv_wj= dv(x_w,y_w,x_j,y_j,fv);
    sv_wj= sv(x_w,y_w,z_w,x_j,y_j,z_j,fv);
    Pwj=P_expt(gv,fv,W,H,X,Y,t,es, dv_wj,sv_wj);

    g=gain(eta,inc,inc_r,fov);
    
    [v1,d1]=point_to_vector(x_i,y_i,z_i,x_w,y_w,z_w);
    Nnorm1=norm_vec_trans(alpha_i,beta_i);
    p1=dot_product(v1,Nnorm1);
    [v2,d2]=point_to_vector(x_w,y_w,z_w,x_i,y_i,z_i);
    Nnorm2=norm_vec_receiver(alpha_w,beta_w);
    p2=dot_product(v2,Nnorm2);
    [v3,d3]=point_to_vector(x_w,y_w,z_w,x_j,y_j,z_j);
    Nnorm3=norm_vec_receiver(alpha_w,beta_w);
    p3=dot_product(v3,Nnorm3);
    [v4,d4]=point_to_vector(x_j,y_j,z_j,x_w,y_w,z_w);
    Nnorm4=norm_vec_receiver(alpha_j,beta_j);
    p4=dot_product(v4,Nnorm4);
    
    digits(2);
    dm=((d1+d3)/c);
    dm=vpa(dm);
    dm=double(subs(dm));
    m_HnLoS= abs(((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g*Piw*Pwj)/((d1^2)*(d3^2)*d1*d2*d3*d4));
    m_HnLoS=vpa(m_HnLoS);
    m_HnLoS=double(subs(m_HnLoS));
end



%función para calcular el HLoS%
function [m_HLoS,dm] = HLoS_direct(x_i,y_i,z_i,x_j,y_j,z_j,Ap,eta,alpha_i,alpha_j,beta_i,beta_j,incidencia,incidencia_r,m,fov,gv,fv,W,H,X,Y,t,es,c)
    dv_ij= dv(x_i,y_i,x_j,y_j,fv);
    sv_ij= sv(x_i,y_i,z_i,x_j,y_j,z_j,fv);
    Pij=P_expt(gv,fv,W,H,X,Y,t,es, dv_ij,sv_ij);
    [v1,d1]=point_to_vector(x_i,y_i,z_i,x_j,y_j,z_j);
    Nnorm1=norm_vec_trans(alpha_i,beta_i);
    p1=dot_product(v1,Nnorm1);
    [v2,d2]=point_to_vector(x_j,y_j,z_j,x_i,y_i,z_i);
    Nnorm2=norm_vec_receiver(alpha_j,beta_j);
    p2=dot_product(v2,Nnorm2);
    g=gain(eta,incidencia,incidencia_r,fov);
    digits(2);
    dm= d1/c;
    dm=vpa(dm);
    dm=double(subs(dm));

    if (incidencia>=0) && (incidencia<=2*fov)
        m_HLoS=abs(((m+1)*Ap/(2*3.1416*d1^2))*(p1^m/d1)*(p2/d2)*g* Pij);
        m_HLoS=vpa(m_HLoS);
        m_HLoS=double(subs(m_HLoS));
    else
        m_HLoS=0;
    end
end


% Function to calculate the  Hscat for 40 scatering point
function [m_total_Hscat,dm_total] = H_scater(x_i,y_i,z_i,x_j,y_j,z_j,Ap,m,f,g,gymma,kr,km,ks,p,N,theta_ij,c,alpha_i,beta_i)
    [v1,dij]=point_to_vector(x_i,y_i,z_i,x_j,y_j,z_j);
    dm_total=zeros(1,(N+4));
    m_total_Hscat=zeros(1,N);

    for i =1:N
        Rr=0.5;
        rn=Rr*rand(1,1);
        theta_sn_j=randi([-180 180]);
        B_i_sn=Bisn(theta_sn_j,theta_ij);

        xs=rn*cosd(B_i_sn);
        ys=rn*sind(B_i_sn);
        zs=rn*cosd(theta_sn_j);
        phi_i_sn_radian=phi_scater(xs,ys,zs,x_i,y_i,z_i,alpha_i,beta_i);
        phi_i_sn=rad2deg(phi_i_sn_radian);

        di_sn=sqrt(rn^2+dij^2-2*rn*dij*cosd(B_i_sn));
        Di_j=di_sn+rn;
        Gn=Gain_n(f,g,gymma,phi_i_sn,kr,km,ks,p,N);
            
        digits(2);
        dm=(Di_j/c);
        dm=vpa(dm);
        dm=double(subs(dm));
        dm_total(i)=+dm;

        if (theta_sn_j>=-180) && (theta_sn_j<=180)
            Hscat=abs(((m+1)*Ap*Gn/(2*3.1416*Di_j^2))*(cosd(phi_i_sn))^m*cosd(theta_sn_j));
            Hscat=vpa(Hscat);
            Hscat=double(subs(Hscat));
        else
            Hscat=0;
        end
        m_total_Hscat(i)=+Hscat;
    end
end


%Function to calculate the  Pij for shadowing model
function Pij = P_expt(gv,fv,W,H,X,Y,t,es,d_v,s_v)
    syms w h x y p E

    if(gv(1)>=2*d_v) &&(gv(2)>=s_v)
        w_int=int(gv(1),w,0,W);
        h_int=int(gv(2),h,0,H);
        A=[w_int h_int];
        x_int=int(fv(1),x,0,X);
        y_int=int(fv(2),y,0,Y);
        B=[x_int y_int];
        exp_value=dot(A,B);
        f=p*t;
        est=-es*exp_value;
        d=vpa(subs(f,p,est),8);
        f=exp(E);
        Pij=vpa(subs(f,E,d),4);
    else
        Pij=0;
    end
end


%function to calculate the d(xv,yv)  
function d_v = dv(x1,y1,x2,y2,fv)
    d_v= abs((y1-y2)*fv(1)-(x1-x2)*fv(2)-x2*y1+x1*y2)/sqrt((y1-y2)^2+(x1-x2)^2);
end


%function to calculate the s(xv,yv) 
function s_v = sv(x1,y1,z1,x2,y2,z2,fv)
    if(z1<=z2)
        s_v = ((y1-y2)^2+(x1-x2)^2+(fv(1)-x1)^2+(fv(2)-y1)^2-[(fv(1)-x2)^2+(fv(2)-y2)^2]/2*sqrt((y1-y2)^2+(x1-x2)^2))+z1;
    else
        s_v = ((y1-y2)^2+(x2-x1)^2+(fv(1)-x2)^2+(fv(2)-y2)^2-[(fv(1)-x1)^2+(fv(2)-y1)^2]/2*sqrt((y1-y2)^2+(x2-x1)^2))+z2;
    end
end


%funcion para calcular la ganancia
function g = gain(eta,incide,incide_r,fov)
    if (incide_r<=2*fov) && (incide_r>=0)
        g = (eta^2)/(sin(incide)^2);
    else
        g = 0;
    end
end


% funcion para calcular el producto punto
function f = dot_product(a,b)
     f=dot(a,b);
end

% Calculo del vector normal para receptores
function m = norm_vec_receiver(alpha,beta)
     m = [cosd(alpha)*sind(beta),sind(alpha)*sind(beta),cosd(beta)];
end


% Calculo del vector normal para transmisores
function n = norm_vec_trans(alpha,beta)
     n = [cosd(alpha)*sind(beta),sind(alpha)*sind(beta),-cosd(beta)];
end


%función para calcular vector y distancia entre dos puntos
function [Vec, lenght] = point_to_vector(x1,y1,z1,x2,y2,z2)
     Vec = [x2-x1,y2-y1,z2-z1];
     lenght = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
end


%calculation of Gn gain
function Gn = Gain_n(f,g,gymma,phi_i_sn,kr,km,ks,p,N)
    p_mie = pmie(f,g,phi_i_sn);
    p_ray = pray(gymma,phi_i_sn);
    p_total = (kr/ks)*p_ray+(km/ks)*p_mie;
    f_scat = p_total*sind(phi_i_sn);
    Gn = p*f_scat/N;
end


%calculation of pmie
function p_mie = pmie(f,g,phi_i_sn)
    p_mie = (1-g^2/4*pi)*(1/(1+g^2-2*g*cosd(phi_i_sn))^1.5+ f*3*(cosd(phi_i_sn))^2-1/2*(1+g^2)^1.5);
end


%calculation of pray
function p_ray = pray(gymma,phi_i_sn)
    p_ray = 3*[1+3*gymma+(1-gymma)*(cosd(phi_i_sn))^2]/(16*pi*(1+2*gymma));
end


%scatering angles phi and theta
function phi_i_sn = phi_scater(x1,y1,z1,x2,y2,z2,alpha,beta)
    v = [x1-x2,y1-y2,z1-z2];
    Ntilt = [cosd(alpha)*sind(beta),sind(alpha)*sind(beta),-cosd(beta)];
    d_p = dot_product(v,Ntilt);
    d = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    phi_i_sn = acos(d_p/d);
end


%scatering angle beta 
function B_i_sn = Bisn(theta_sn_j,theta_ij)
    if(theta_ij<theta_sn_j)
        B_i_sn = theta_sn_j-theta_ij;
    else
        B_i_sn = theta_ij-theta_sn_j;
    end
end


% function to plot for one impulse
function plott = plottt(y1,t,dij,D)
    d = y1;
    idx = d == Inf; % find Inf
    d(idx) = dij;
    plott = plot3(t,D*ones(size(t)),d);
end


% Distance function
function D = DIST(x1,y1,z1,x2,y2,z2)
     D = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
end
   