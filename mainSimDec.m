

%%Initialize
clear
clc
close all

tic

%%Define Key Locations
% Area of Interest (AOI)
AOIs.center_lat = 35.4212; % Siachen Glacier latitude
AOIs.center_long = 77.1090; % Siachen Glacier longitude
AOIs.radius_km = 200; % Radius of AOI in km

% Area of Interest (AOI)
AOIb.center_lat = 35.7128; % Blator Glacier latitude
AOIb.center_long = 76.5133; % Blatori Glacier longitude
AOIb.radius_km = 200; % Radius of AOI in km






% Ground Station (GS)
GS.lat = 33.6844; % Islamabad latitude
GS.long = 73.0479; % Islamabad longitude


%Mission Moded
% Mainttain mission modes 
% 0---- Detumble
% 1---- Active
% 2-----SafeMode

Mission_mode= 0; % 



% EPS Parameters
% EPS Parameters
solar_panel_area = 0.1; % Solar panel area (m^2)
panel_efficiency = 0.28; % Solar panel efficiency
solar_constant = 1361; % Solar constant (W/m^2)
battery_capacity = 50; % Battery capacity (Watt-hours)
battery_storage = battery_capacity * 3600; % Convert to Joules
initial_battery_level = battery_storage; % Initial battery charge (Joules)
safe_mode_threshold = 0.1 * battery_storage; % 10% battery level threshold for safe mode
% Subsystem Power Consumption (W)
power_baseline = 5; % Satellite bus (OBC, sensors)
power_comm = 2; % Communication system
power_payload = 6; % Payload
power_adcs = 4; % ADCS during detumbling
power_safe_mode = 3; % Safe mode power consumption

power_total=power_baseline ;


%Communication  Parameters
data_to_transmit= 0; % payload data
data_remaining = 0;  % Remaining data to transmit initalise after image capture
packet_idx=0;
data_transmitted=0;
%downlink_bandwidth=2*1e6; % Bandwidth 
comm_time_step = 10;
data_rate = 9600; % Data rate (bits per second)
packet_size = 256; % Packet size (bytes)
preamble = uint8([126, 126]); % AX.25 preamble flag (0x7E)
% AX.25 Header Fields
source_address = 'SAT001'; % 6-character satellite ID
destination_address = 'GND001'; % 6-character ground station ID
control_field = uint8(0x03); % Control field (default AX.25 value)
protocol_id = uint8(0xF0); % Protocol ID for unconnected data transfer




%Payload
Payload.imageData = {}; % Stores images
Payload.Sent = {}; % Status of whether the image has been sent
 % 
 % Payload.imageData = {};
 % Payload.Sent = 0; % image is delivered or not
%% Initialize Log
comm_status = struct('time', [], 'data_transmitted', [], 'data_remaining', [], 'downlink_times', []);
downlink_start = NaN; % To track start of downlink




%%%Globals
global BI BB m I Is invI mu lastMagUpdate nextMagUpdate lastSensorUpdate 
global nextSensorUpdate BfieldMeasured pqrMeasured ptpMeasured BfieldNav pqrNav ptpNav
global BfieldNavPrev pqrNavPrev ptpNavPrev current Ir1Bcg Ir2Bcg Ir3Bcg n1 n2 n3
global maxSpeed maxAlpha Ir1B Ir2B Ir3B rwalphas
global fsensor MagFieldBias AngFieldBias EulerBias R Amax lmax CD
global MagFieldNoise AngFieldNoise EulerNoise IrR Jinv
 
%%%%Initialize Nav Filter
BfieldNavPrev = [-99;0;0];
pqrNavPrev = [0;0;0];
ptpNavPrev = [0;0;0];

%%%%Simulation of a Low Earth Satellite
disp('Simulation Started')


nextMagUpdate = 1;
lastMagUpdate = 0;

%%%Get Planet Parameters
planet

%%%Get mass and inertia properties
inertia

%%%Initial Conditions Position and Velocity
altitude = 600*1000; %%meters
x0 = R + altitude;
y0 = 0;
z0 = 0; 
xdot0 = 0;
inclination =35*pi/180;

semi_major = norm([x0;y0;z0]);
vcircular = sqrt(mu/semi_major);
ydot0 = vcircular*cos(inclination);
zdot0 = vcircular*sin(inclination);

%%%Intitial Conditions for Attitude and Angular Velocity
phi0 = 0;
theta0 = 0;
psi0 = 0;
ptp0 = [phi0;theta0;psi0];
q0123_0 = EulerAngles2Quaternions(ptp0);

%Initital  Angular Velocity

% p0 = 2;
% q0 = 2;
% r0 = 2;



p0 = 1.8;
q0 = -0.5;
r0 = 0.3;
%%%Initial conditions of my reaction wheels
% w10 = 0;
% w20 = 0;
% w30 = 0;



state = [x0;y0;z0;xdot0;ydot0;zdot0;q0123_0;p0;q0;r0];

%%%Need time window
period = 3*pi/sqrt(mu)*semi_major^(3/2);
number_of_orbits = 1;
tfinal = period*number_of_orbits;
%tfinal = 100;
next = 10;
timestep = 1.0;
tout = 0:timestep:tfinal;
stateout = zeros(length(tout),length(state));


%%%%%%%%

% Initialize variables EPS
battery_level = zeros(length(tout), 1); % Battery level (Joules)
solar_power = zeros(length(tout), 1); % Solar power output (W)
energy_consumed = zeros(length(tout), 1); % Energy consumed (Joules)
subsystem_mode = cell(length(tout), 1); % Subsystem operation mode
in_sunlight = true(length(tout), 1); % Sunlight or eclipse state
battery_level(1) = initial_battery_level; % Initial battery level
is_safe_mode = false; % Safe mode flag





%%%This is where we integrate the equations of motion

%%%Loop through time to integrate
BxBout = 0*stateout(:,1);
ByBout = BxBout;
BzBout = BxBout;
BxBm = 0*stateout(:,1);
ByBm = BxBout;
BzBm = BxBout;
pqrm = zeros(length(tout),3);

ptpm = zeros(length(tout),3);
ptpN = 0*ptpm;

BxBN = 0*stateout(:,1);
ByBN = BxBout;
BzBN = BxBout;
pqrN = zeros(length(tout),3);

ix = 0*stateout(:,1);
iy = ix;
iz = ix;

%rwa = 0*ptpm;

%%%Sensor Parameters
lastSensorUpdate = 0;
sensor_params

%%%%Call the Derivatives Routine to initialize vars
k1 = Satellite(tout(1),state);

%%%Print Next
lastPrint = 0;
for idx = 1:length(tout)
    %%%Save the current state
    stateout(idx,:) = state';
    
    %%%Save the Current
    ix(idx) = current(1);
    iy(idx) = current(2);
    iz(idx) = current(3);
    
    % %%%%Save reaction wheel acceleration
    % rwa(idx,:) = rwalphas';
    
    %%%Save the magnetic field
    BxBout(idx) = BB(1);
    ByBout(idx) = BB(2);
    BzBout(idx) = BB(3);


%%%% tracing of intermediat results

 Btemp(idx,1) = BI(1);
    Btemp(idx,2) = BI(2);
    Btemp(idx,3) = BI(3);



    %%%Save the Measured Magnetic field
    BxBm(idx) = BfieldMeasured(1);
    ByBm(idx) = BfieldMeasured(2);
    BzBm(idx) = BfieldMeasured(3);
    %%%Save the Nav magnetic field
    BxBN(idx) = BfieldNav(1);
    ByBN(idx) = BfieldNav(2);
    BzBN(idx) = BfieldNav(3);
    %%%THe actual pqr truth signal is embedded in
    %%%The state vector.
    %%%Save the measured pqr signal
    pqrm(idx,:) = pqrMeasured';
    %%%Save the Nav pqr signal
    pqrN(idx,:) = pqrNav';
    %%%Save ptp
    ptpm(idx,:) = ptpMeasured';
    ptpN(idx,:) = ptpNav';


    %%%%% Solar Panels
    % Compute current orbital position (simplified circular orbit)
    true_anomaly = mod((tout(idx) / period) * 2 * pi, 2 * pi); % True anomaly (radians)

    % Determine sunlight or eclipse (simplified shadow model)
    if true_anomaly > pi % In Earth's shadow
        in_sunlight(idx) = false;
        solar_power(idx) = 0; % No solar power in eclipse
    else
        in_sunlight(idx) = true;
        solar_power(idx) = solar_constant * solar_panel_area * panel_efficiency; % Solar power (W)
    end


    
   %%Detumbling if angular velocity is high wait till detumble 
   
    if (max(abs(Quaternions2EulerAngles(state(7:10)'))) > 0.001) && ( is_safe_mode == false)
        Mission_mode= 0;
        power_total = power_baseline + power_adcs;
        subsystem_mode{idx} = 'Detumbling';
        fprintf("Detumbling \n")
    else
         Mission_mode= 1;% Active Mission Mode
     %    fprintf("Active \n")
    end

     if is_safe_mode
        power_total = power_safe_mode;
        subsystem_mode{idx} = 'Safe Mode';
        Mission_mode=2;
     end

if( Mission_mode == 1)

    
    % ------------Communication

        [lat,long]=eci_to_latlong(state(1,1),state(2,1),state(3,1));
        % fprintf('lat: %.1f long%.1f.\n', lat, long);
          distance_to_gs = haversine_distance(lat, long, GS.lat, GS.long);
            if distance_to_gs <= 250 % Assuming 250 km communication range
                 fprintf('Time: %.1f s - CubeSat is over the Ground Station.\n', idx);
                 power_total = power_baseline + power_comm; %Powers Consumption
                  subsystem_mode{idx} = 'Communication';
                  for i = 1:50 % Add a rate simluatio nstep
                    if( data_remaining > 0)
        
        
                        % Extract payload for the current packet
                        start_idx = (packet_idx - 1) * packet_size + 1;
                        end_idx = min(packet_idx * packet_size, length(image_vector));
                        payload_data = image_vector(start_idx:end_idx);
        
                        % Build AX.25 Frame
                        lat_long = typecast(single([lat, long]), 'uint8'); % Convert lat/long to bytes
                        header = buildAX25Header(source_address, destination_address, control_field, protocol_id);
                        frame = [preamble, header, lat_long, payload_data'];
        
                        % Add Frame Check Sequence (FCS)
                        fcs = calculateCRC(frame);
                        frame = [frame, fcs];
        
                       packet_idx = packet_idx+1;
                       data_remaining = data_remaining -1;
                       if (data_remaining < 1)
                           length(Payload.imageData)
                           Payload.Sent{1}=true;
                           %%%%%Keep unsent data only
                            toKeep = ~[Payload.Sent{:}];
                            Payload.imageData = Payload.imageData(toKeep);
                            Payload.Sent = Payload.Sent(toKeep);
                            length(Payload.imageData);
                       end
                            data_transmitted=data_transmitted+1;

                    % Log Communication Status
                    comm_status.time = [comm_status.time; idx];
                    comm_status.data_transmitted = [comm_status.data_transmitted; data_transmitted];
                    comm_status.data_remaining = [comm_status.data_remaining; data_remaining];
        
                     fprintf('Time: %.1f s - Data Pakect Transmitted: %.2f MB, Data Remaining: %.2f MB\n',idx, data_transmitted,  data_remaining);
        
                    end
        

                   end
        
               end
        
 
%----------Payload        
                   
            distance_to_gs = haversine_distance(lat, long, AOIb.center_lat, AOIb.center_long);
            if distance_to_gs <= AOIb.radius_km % Assuming 250 km communication range
                fprintf('Time: %.1f s - CubeSat is over Baltoro.\n', idx);
                    fprintf('Imaging Mode: Capturing image...\n');
                    image = struct('timestamp', datetime, 'data', rand(256, 256)); % Simulated image
                    Payload.imageData{end+1} = image;
                    Payload.Sent{end+1}= false;

                    power_total = power_baseline + power_payload; % Payload Power
                    subsystem_mode{idx} = 'Payload Operation';
            end 
        
        
           distance_to_gs = haversine_distance(lat, long, AOIs.center_lat, AOIs.center_long);
            if distance_to_gs <= AOIs.radius_km % Assuming 250 km communication range
                fprintf('Time: %.1f s - CubeSat is over Siachen.\n', idx);
                    fprintf('Imaging Mode: Capturing image...\n');
                    image = struct('timestamp', datetime, 'data', rand(256, 256)); % Simulated image
                    Payload.imageData{end+1} = image;
                    Payload.Sent{end+1}= false;

                     power_total = power_baseline + power_payload; % Payload Power
                    subsystem_mode{idx} = 'Payload Operation';
            end 
        
        %%%% Convert data in to frames
            if((length(Payload.imageData)>=1) && (data_remaining <= 0))
                image_vector = reshape(Payload.imageData{1}.data, [], 1); % Flatten the image into a 1D array
                num_packets = ceil(length(image_vector) / packet_size); % Number of packets
                data_remaining= num_packets; % remaining data to sent
                packet_idx=1; %Packet to be sent
            end       

end
%----------EPS-----------------
  % Energy consumption
    energy_consumed(idx) = power_total * timestep;

    % Update battery level
    if idx > 1
        if in_sunlight(idx)
            % Battery charges in sunlight
            net_power = solar_power(idx) - power_total;
            battery_level(idx) = min(battery_level(idx-1) + net_power * timestep, battery_storage);
        else
            % Battery discharges in eclipse
            battery_level(idx) = max(battery_level(idx-1) - power_total * timestep, 0);
        end
    end

    % Check for safe mode activation
    if battery_level(idx) < safe_mode_threshold
        is_safe_mode = true;
    else
         is_safe_mode = false;
    end








    % 
    %%%%Then we make our 4 function calls for the RK4
    k1 = Satellite(tout(idx),state);
    k2 = Satellite(tout(idx)+timestep/2,state+k1*timestep/2);
    k3 = Satellite(tout(idx)+timestep/2,state+k2*timestep/2);
    k4 = Satellite(tout(idx)+timestep,state+k3*timestep);
    k = (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    state = state + k*timestep;
    
    
end
%%%Save original State
stateout_original = stateout;
disp('Simulation Complete')
toc

%%Convert state to kilometers
stateout(:,1:6) = stateout_original(:,1:6)/1000;

%%%Extract the state vector
xout = stateout(:,1);
yout = stateout(:,2);
zout = stateout(:,3);
q0123out = stateout(:,7:10);
ptpout = Quaternions2EulerAngles(q0123out);
pqrout = stateout(:,11:13);
%w123 = stateout(:,14:16);


%%%Make an Earth
[x_sphere, y_sphere, z_sphere] = sphere(50);

%%%Plot 3D orbit
fig = figure();
set(fig,'color','white')

plot3(xout,yout,zout,'r','LineWidth',4)
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
hold on
% Plot Earth as a Sphere
surf(1e-3*R * x_sphere, 1e-3*R * y_sphere, 1e-3*R * z_sphere, ...
    'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
axis equal


% Plot AOI as a Point
[aoi_x, aoi_y, aoi_z] = latlong_to_cartesian(AOIs.center_lat, AOIs.center_long, 1e-3*R);
scatter3(aoi_x, aoi_y, aoi_z, 100, 'g', 'filled', 'DisplayName', 'Area of Interest');
hold on;

[aoi_x, aoi_y, aoi_z] = latlong_to_cartesian(AOIb.center_lat, AOIb.center_long, 1e-3*R);
scatter3(aoi_x, aoi_y, aoi_z, 100, 'g', 'filled', 'DisplayName', 'Area of Interest');
hold on;

% Plot Ground Station
[gs_x, gs_y, gs_z] = latlong_to_cartesian(GS.lat, GS.long, 1e-3*R);
scatter3(gs_x, gs_y, gs_z, 100, 'y', 'filled', 'DisplayName', 'Ground Station');

%view([-227 24]);
 hold off;

view([-186 31])



%%%plot Euler Angles
fig4 = figure();
set(fig4,'color','white')
p1 = plot(tout,ptpout*180/pi,'-','LineWidth',2);
hold on
p2 = plot(tout,ptpm*180/pi,'-s','LineWidth',2);
p3 = plot(tout,ptpN*180/pi,'--','LineWidth',2);
grid on
xlabel('Time (sec)')
ylabel('Euler Angles (deg)')
legend('Phi','Theta','Psi')
legend([p1(1),p2(1),p3(1)],'Actual','Measured','Nav')





% Convert battery level to Watt-hours for easier interpretation
battery_level_wh = battery_level / 3600;

% Visualization
figure;

% Plot Battery Level
subplot(4, 1, 1);
plot(tout / 3600, battery_level_wh, 'b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Battery Level (Wh)');
title('Battery Level Over Time');
grid on;

% Plot Solar Power
subplot(4, 1, 2);
plot(tout / 3600, solar_power, 'r', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Solar Power (W)');
title('Solar Power Output Over Time');
grid on;

% Plot Energy Consumption
subplot(4, 1, 3);
plot(tout / 3600, energy_consumed / timestep, 'g', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Energy Consumed (W)');
title('Energy Consumption Over Time');
grid on;

% Plot Subsystem Mode
subplot(4, 1, 4);
hold on;
for idx = 1:length(tout)
    if strcmp(subsystem_mode{idx}, 'Detumbling')
        plot(tout(idx) / 3600, 1, 'ro');
    elseif strcmp(subsystem_mode{idx}, 'Communication')
        plot(tout(idx) / 3600, 2, 'bo');
    elseif strcmp(subsystem_mode{idx}, 'Payload Operation')
        plot(tout(idx) / 3600, 3, 'go');
    else
        plot(tout(idx) / 3600, 0, 'ko');
    end
end
xlabel('Time (hours)');
ylabel('Subsystem Mode');
title('Subsystem Mode Over Time');
yticks([0, 1, 2, 3]);
yticklabels({'Safe Mode', 'Detumbling', 'Communication', 'Payload Operation'});
grid on;


toc


function [x, y, z] = latlong_to_cartesian(lat, long, radius)
    % Convert latitude and longitude to Cartesian coordinates
    x = radius * cosd(lat) * cosd(long);
    y = radius * cosd(lat) * sind(long);
    z = radius * sind(lat);
end
function distance = haversine_distance(lat1, long1, lat2, long2)
    % Calculate the great-circle distance between two points on a sphere
    R = 6371; % Earth's radius in km
    dlat = deg2rad(lat2 - lat1);
    dlong = deg2rad(long2 - long1);
    a = sin(dlat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlong/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distance = R * c;
end
function [lat, long] = eci_to_latlong(x, y, z)
    % Convert ECI coordinates to latitude and longitude
   
    rho = norm([x;y;z]);

    thetaE = acos(z/rho);
    psiE = atan2(y,x);

    
    lat = 90 - rad2deg(thetaE); %90-thetaE*180/p;
    long = rad2deg(psiE);%*180/pi;




    % lat = asind(z / rho);
    % long = atan2d(y, x);
   
end
function header = buildAX25Header(src, dest, control, pid)
    % AX.25 Header Construction
    header = [encodeCallsign(dest), encodeCallsign(src), control, pid];
end

function encoded = encodeCallsign(callsign)
    % Encode a 6-character callsign into AX.25 format
    callsign = pad(callsign, 6, 'right', ' '); % Ensure callsign is 6 characters
    encoded = uint8(callsign) * 2; % Shift ASCII values left by 1 bit
end

function crc = calculateCRC(data)
    % CRC-16-CCITT Calculation (FCS)
    polynomial = uint16(0x1021);
    crc = uint16(0xFFFF);
    for byte = data
        crc = bitxor(crc, bitshift(uint16(byte), 8));
        for bit = 1:8
            if bitand(crc, uint16(0x8000))
                crc = bitxor(bitshift(crc, 1), polynomial);
            else
                crc = bitshift(crc, 1);
            end
        end
    end
    crc = uint8([bitshift(crc, -8), bitand(crc, uint16(0xFF))]); % Convert to 2 bytes
end


