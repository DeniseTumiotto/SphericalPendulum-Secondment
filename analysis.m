% % [sols, params] = readAll();
% [sols, params] = readUpto(2);
% % evaluate error
% 
% % error = evalErr(sols, params);
% % plotErr(error)
% 
% for i = 0:floor(size(sols,1)/2)
%     my_dist1 = riemannianDistance(sols(2*i+1:2*i+2), params(2*i+1:2*i+2));
%     
%     figure()
%     time = linspace(params{2*i+1}.t0, params{2*i+1}.T, params{2*i+1}.N_TIME);
%     plot(time, my_dist1, 'b-', 'LineWidth', 3)
%     title('Riemaniann Distance', 'FontSize', 16)
%     xlabel('time in s', 'FontSize', 16)
%     ylabel('distance in m', 'FontSize', 16)
%     grid on
% 
% end
% 
% % animationSphere(sols(1), params(1), 100000);


%% linear analysis

% [s, p] = readUpto(6);
% 
% my_dist_dmp = euclidDistance(s{1}, s{2});
% my_dist = euclidDistance(s{3}, s{4});
% 
% my_dist_lie = riemannianDistance(s(5:6), p(5:6));
% 
% figure()
% % time = linspace(p{1}.t0, p{1}.T, p{1}.N_TIME);
% time = linspace(p{1}.t0, p{1}.te, p{1}.nTime);
% plot(time, my_dist, 'b-', 'LineWidth', 3)
% title('Euclidean Distance', 'FontSize', 16)
% xlabel('time in s', 'FontSize', 16)
% ylabel('distance in m', 'FontSize', 16)
% grid on
% 
% figure()
% plot(time, my_dist_dmp, 'b-', 'LineWidth', 3)
% title('Euclidean Distance', 'FontSize', 16)
% xlabel('time in s', 'FontSize', 16)
% ylabel('distance in m', 'FontSize', 16)
% grid on
% 
% figure()
% time = linspace(p{5}.t0, p{5}.T, p{5}.N_TIME);
% plot(time, my_dist_lie, 'r-', 'LineWidth', 3)
% title('Euclidean Distance', 'FontSize', 16)
% xlabel('time in s', 'FontSize', 16)
% ylabel('distance in m', 'FontSize', 16)
% grid on

%% Problem eigenvalues

clearvars
clc

[s, p] = readUpto(2);

% num = @(x, y) x*transpose(y);
% denom = @(x) 1+transpose(x)*x;
% dhdy = @(x, y, crrD) [num(x,x)+num(y,y)/crrD num(y,x)/crrD;
%         num(x,y)/crrD num(x,x)/crrD];
% 
for k = 1:2
q = s{k}(1:3,:);
w = s{k}(4:6,:);
sph = vec2sph(q, w);
% g = 9.81;
% d = p{k}.damp;
% e3 = [0; 0; 1];
% 
% dt = p{k}.dt;
% for i = 1:size(q,2)-1
%     y = [q(:,i); w(:,i)];
%     A = [skw(w(:,i)) zeros(3); g*skw(e3) -d*eye(3)];
%     Anew = [skw(w(:,i+1)) zeros(3); g*skw(e3) -d*eye(3)];
%     crrDenom = denom(w(:,i));
%     aux = dhdy(q(:,i),w(:,i),crrDenom);
%     P = eye(6) - aux;
%     DP = ((eye(6)-dhdy(q(:,i+1),w(:,i+1),denom(w(:,i+1))))-P)/dt;
%     DA = (Anew - A)/dt;
%     a = DP * A * y;
%     b = P * DA * y;
%     c = P * A;
% 
%     DF = a + b + c;
%     DF_symm = 0.5 * (DF + transpose(DF));
% 
%     lambda = eig(DF_symm);
% 
%     disp(lambda)
% 
%     if any(lambda > 0)
%         index_pos = find(lambda > 0);
%         if any(lambda(index_pos) < 1e-15)
%             index_small = find(lambda(index_pos) < 1e-15);
%             if size(index_small,1) ~= size(index_pos,1)
%                 disp('positive eigenvalue!')
%             end
%         else
%             disp('positive eigenvalue!')
%         end
%     end
% end
end












