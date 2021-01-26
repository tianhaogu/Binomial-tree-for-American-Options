%% Read data file & extract parameters & add necessary variables & Convergence test
apple_data = csvread('AAPL_Option_Data.csv');
interest_data = csvread('Interest_rates.csv');
volatility_data = csvread('HV.csv');
num_contract = size(apple_data,1);
new = ones(num_contract,1); diff = zeros(num_contract,1);
apple_data = [apple_data,new]; 
exercise_max_col = zeros(num_contract,1); binomial_max_col = zeros(num_contract,1); exercise_time_col = zeros(num_contract,1); is_plot = 0;
for i = 1:1:num_contract
   strike_price = apple_data(i,5);
   ttm = apple_data(i,2);
   num_step = 100;
   judge1 = 0;
   volatility = 0;
   for j = 2:1:9
       if ttm <= volatility_data(1,j) && judge1 == 0
           volatility = volatility_data(2,j);
           judge1 = judge1 + 1;
       end
       if j == 9 && judge1 == 0
           volatility = volatility_data(2,j);
       end
   end
   stock_price = apple_data(i,6);
   judge2 = 0;
   interest = 0;
   for k = 16:(-1):1
       if ttm <= interest_data(k,1) && judge2 == 0
           interest = interest_data(k,2)/100;
           judge2 = judge2 + 1;
       end
       if k == 1 && judge2 == 0
           interest = interest_data(k,2)/100;
       end
   end
   option_price = 0; exercise_max = 0; binomial_max = 0; exercise_time = 0;
   if apple_data(i,3) == 1 && apple_data(i,4) == 1
       [option_price,exercise_max,binomial_max, exercise_time] = American_put(strike_price, ttm, num_step, ...
                                                                              volatility, stock_price, interest, is_plot);
       apple_data(i,end) = option_price;
       diff(i,1) = (apple_data(i,end) - apple_data(i,9)) / (apple_data(i,9));
       exercise_max_col(i,1) = exercise_max;
       binomial_max_col(i,1) = binomial_max;
       exercise_time_col(i,1) = exercise_time;
   end
   if apple_data(i,3) == 0 && apple_data(i,4) == 1
       [option_price,exercise_max,binomial_max, exercise_time] = American_call(strike_price, ttm, num_step, ...
                                                                               volatility, stock_price, interest);
       apple_data(i,end) = option_price;
       diff(i,1) = (apple_data(i,end) - apple_data(i,9)) / (apple_data(i,9));
       exercise_max_col(i,1) = exercise_max;
       binomial_max_col(i,1) = binomial_max;
       exercise_time_col(i,1) = exercise_time;
   end
   val_col = zeros(150,1);
   if i == 2 || i == 6
       [val_col] = Convergence_test(strike_price, ttm, volatility, stock_price, interest, apple_data(i,3));
       figure;
       x = apple_data(1:150,1); y = val_col;
       plot(x,y,'b');  set(gca,'fontsize',14);
       leg = legend('calculated price','location','best'); set(leg,'FontSize',18);
       xlabel('number of steps','fontsize',18); ylabel('option price','fontsize',18);
       title('Convergence test---Calculated option price using incremental number of steps','fontsize',20);
   end
end 
apple_data = [apple_data,diff];
apple_data = [apple_data,exercise_max_col];
apple_data = [apple_data,binomial_max_col];
apple_data = [apple_data,exercise_time_col];

%% Plot of comparison between calculated and actual option price
figure;
x = apple_data(1:num_contract,1);
y1 = apple_data(1:num_contract,24);
y2 = apple_data(1:num_contract,9);
scatter(x,y1,'b');
hold on;
plot(x,y2,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('option price','fontsize',18);
title('Option price comparison','fontsize',20);
xlim([1,num_contract]); ylim([146,154]);

%% Plot of accuracy
figure;
y3 = apple_data(1:num_contract,25);
scatter(x,y3,'g'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('error','fontsize',18);
title('Error between actual price and calculated price','fontsize',20);
xlim([1,num_contract]); ylim([-0.03,0.01]);

%% Plot of comparison between calculated call and put option price
apple_put = apple_data(apple_data(:,3)==1,:);
apple_call = apple_data(apple_data(:,3)==0,:);
num_contract_put = size(apple_put,1);
num_contract_call = size(apple_call,1);
x_put = [1:1:num_contract_put];
y_put = apple_put(1:num_contract_put,24);
y_put2 = apple_put(1:num_contract_put,9);
x_call = [1:1:num_contract_call];
y_call = apple_call(1:num_contract_call,24);
y_call2 = apple_call(1:num_contract_call,9);
figure;
subplot(2,1,1);
scatter(x_put,y_put,'b');
hold on;
plot(x_put,y_put2,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('put option price','fontsize',18);
title('Put option price comparison','fontsize',20);
xlim([1,num_contract_put]); ylim([146,154]);
subplot(2,1,2);
scatter(x_call,y_call,'b');
hold on;
plot(x_call,y_call2,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('call option price','fontsize',18);
title('Call option price comparison','fontsize',20);
xlim([1,num_contract_call]); ylim([146,154]);

%% Plot of number/proportion of early exercise of both call and put options
y_call_exercise = apple_call(1:num_contract_call,26);
y_put_exercise = apple_put(1:num_contract_put,26);
figure;
subplot(2,1,1);
scatter(x_call,y_call_exercise,'b'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('Number of early exercise','fontsize',18);
title('Call number of early exercise','fontsize',20);
xlim([1,num_contract_call]); ylim([0,500]);
subplot(2,1,2);
scatter(x_put,y_put_exercise,'b'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('Number of early exercise','fontsize',18);
title('Put number of early exercise','fontsize',20);
xlim([1,num_contract_put]); ylim([0,2500]);

%% Plot of comparisonb among different time to maturity (split into four intervals)
apple_short = apple_data(apple_data(:,2)<=30,:);
apple_medium1 = apple_data(apple_data(:,2)>30 & apple_data(:,2)<=60,:);
apple_medium2 = apple_data(apple_data(:,2)>60 & apple_data(:,2)<=185,:);
apple_long = apple_data(apple_data(:,2)>185,:);
num_short = size(apple_short,1);
num_medium1 = size(apple_medium1,1);
num_medium2 = size(apple_medium2,1);
num_long = size(apple_long,1);
x_short = [1:1:num_short];
x_medium1 = [1:1:num_medium1];
x_medium2 = [1:1:num_medium2];
x_long = [1:1:num_long];
y_short_cal = apple_short(1:num_short,24);
y_short_act = apple_short(1:num_short,9);
y_medium1_cal = apple_medium1(1:num_medium1,24);
y_medium1_act = apple_medium1(1:num_medium1,9);
y_medium2_cal = apple_medium2(1:num_medium2,24);
y_medium2_act = apple_medium2(1:num_medium2,9);
y_long_cal = apple_long(1:num_long,24);
y_long_act = apple_long(1:num_long,9);
figure;
subplot(2,2,1);
scatter(x_short,y_short_cal,'b');
hold on;
plot(x_short,y_short_act,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('option price','fontsize',18);
title('Option price of maturity shorter than 1 month','fontsize',20);
xlim([1,num_short]); ylim([146,154]);
subplot(2,2,2);
scatter(x_medium1,y_medium1_cal,'b');
hold on;
plot(x_medium1,y_medium1_act,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('option price','fontsize',18);
title('Option price of maturity from 1-2 months','fontsize',20);
xlim([1,num_medium1]); ylim([146,154]);
subplot(2,2,3);
scatter(x_medium2,y_medium2_cal,'b');
hold on;
plot(x_medium2,y_medium2_act,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('option price','fontsize',18);
title('Option price of maturity from 2-6 months','fontsize',20);
xlim([1,num_medium2]); ylim([146,154]);
subplot(2,2,4);
scatter(x_long,y_long_cal,'b');
hold on;
plot(x_long,y_long_act,'r'); set(gca,'fontsize',14);
leg = legend('calculated price','actual price','location','best'); set(leg,'FontSize',18);
xlabel('serial number of option contract','fontsize',18); ylabel('option price','fontsize',18);
title('Option price of maturity longer than 6 months','fontsize',20);
xlim([1,num_long]); ylim([146,154]);

%% Plot of time(step) of early exercise, and comparison between call/put option and different time to maturity
apple_short_call = apple_short(apple_short(:,3)==0,:); apple_short_put = apple_short(apple_short(:,3)==1,:);
apple_medium1_call = apple_medium1(apple_medium1(:,3)==0,:); apple_medium1_put = apple_medium1(apple_medium1(:,3)==1,:);
apple_medium2_call = apple_medium2(apple_medium2(:,3)==0,:); apple_medium2_put = apple_medium2(apple_medium2(:,3)==1,:);
apple_long_call = apple_long(apple_long(:,3)==0,:); apple_long_put = apple_long(apple_long(:,3)==1,:);
num_short_call = size(apple_short_call,1); num_short_put = size(apple_short_put,1);
num_medium1_call = size(apple_medium1_call,1); num_medium1_put = size(apple_medium1_put,1);
num_medium2_call = size(apple_medium2_call,1); num_medium2_put = size(apple_medium2_put,1);
num_long_call = size(apple_long_call,1); num_long_put = size(apple_long_put,1);
x_short_call = [1:1:num_short_call]; x_short_put = [1:1:num_short_put];
x_medium1_call = [1:1:num_medium1_call]; x_medium1_put = [1:1:num_medium1_put];
x_medium2_call = [1:1:num_medium2_call]; x_medium2_put = [1:1:num_medium2_put];
x_long_call = [1:1:num_long_call]; x_long_put = [1:1:num_long_put];
y_short_call_exercisetime = apple_short_call(1:num_short_call,28); y_short_put_exercisetime = apple_short_put(1:num_short_put,28);
y_medium1_call_exercisetime = apple_medium1_call(1:num_medium1_call,28); y_medium1_put_exercisetime = apple_medium1_put(1:num_medium1_put,28);
y_medium2_call_exercisetime = apple_medium2_call(1:num_medium2_call,28); y_medium2_put_exercisetime = apple_medium2_put(1:num_medium2_put,28);
y_long_call_exercisetime = apple_long_call(1:num_long_call,28); y_long_put_exercisetime = apple_long_put(1:num_long_put,28);
figure;
subplot(2,2,1);
scatter(x_short_call,y_short_call_exercisetime,'r'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Call Option with maturity less than 1 month','fontsize',20);
xlim([1,num_short_call]); ylim([0,200]);
subplot(2,2,2);
scatter(x_medium1_call,y_medium1_call_exercisetime,'r'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Call Option with maturity from 1-2 months','fontsize',20);
xlim([1,num_medium1_call]); ylim([0,200]);
subplot(2,2,3);
scatter(x_medium2_call,y_medium2_call_exercisetime,'r'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Call Option with maturity from 2-6 months','fontsize',20);
xlim([1,num_medium2_call]); ylim([0,200]);
subplot(2,2,4);
scatter(x_long_call,y_long_call_exercisetime,'r'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Call Option with maturity longer than 6 months','fontsize',20);
xlim([1,num_long_call]); ylim([0,200]);
figure;
subplot(2,2,1);
scatter(x_short_put,y_short_put_exercisetime,'g'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Put Option with maturity less than 1 month','fontsize',20);
xlim([1,num_short_put]); ylim([0,200]);
subplot(2,2,2);
scatter(x_medium1_put,y_medium1_put_exercisetime,'g'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Put Option with maturity from 1-2 months','fontsize',20);
xlim([1,num_medium1_put]); ylim([0,200]);
subplot(2,2,3);
scatter(x_medium2_put,y_medium2_put_exercisetime,'g'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Put Option with maturity from 2-6 months','fontsize',20);
xlim([1,num_medium2_put]); ylim([0,200]);
subplot(2,2,4);
scatter(x_long_put,y_long_put_exercisetime,'g'); set(gca,'fontsize',14);
xlabel('serial number of option contract','fontsize',18); ylabel('early exercise time','fontsize',18);
title('Time of early exercise --- Put Option with maturity longer than 6 months','fontsize',20);
xlim([1,num_long_put]); ylim([0,200]);
hold off;

%% Calculate the proportion of early exercise of both call and put options. Plot binomial tree model 
call_early_exercise = size(apple_call(apple_call(:,26)~=0,:),1) / num_contract_call
put_early_exercise = size(apple_put(apple_put(:,26)~=0,:),1) / num_contract_put
is_plot = 1;
[option_price] = American_put(apple_data(76,5), apple_data(76,2), 40, ...
                              volatility_data(2,5), apple_data(76,6), interest_data(15,2), is_plot);

%% Algorithm of binomial tree model (both put and call) for American option pricing
function [put_value, exercise_max, binomial_max, exercise_time] = American_put(K,T,N,Vol,S,r,is_plot)
  if T > 365
      N = floor(N*1.8);
  end
  dT = T/N/365;
  up = exp(Vol*sqrt(dT));
  down = exp(-Vol*sqrt(dT));
  up_prob = (exp(r*dT)-down) / (up-down);
  down_prob = (up-exp(r*dT)) / (up-down);
  depth = N+1;
  node_price = zeros(depth,1);
  put_price = zeros(depth,1);
  exercise_max = 0; binomial_max = 0; exercise_time = 0; whether_exercise = 0;
  if is_plot == 1
      figure;
      subplot(1,2,1); line([1 depth],[K,K],'color','red');
      hold on;
  end
  for j = depth:(-1):1
      for i = 1:1:j
          node_price(i) = S*up^(j-i)*down^(i-1);
          if j == depth
              put_price(i) = max(K-node_price(i), 0);
          else
              put_price(i) = max((K-node_price(i)), ...
                                 exp(-r*dT)*(up_prob*node_price(i)+down_prob*node_price(i+1)));
              if put_price(i) == K-node_price(i)     % analysis of early exercise
                  exercise_max = exercise_max + 1;
                  if whether_exercise == 0
                      exercise_time = j-1;
                      whether_exercise = whether_exercise+1;
                  end
              else
                  binomial_max = binomial_max + 1;
              end
          end
      end
      whether_exercise = 0;
      if is_plot == 1      % plot of binomial tree model
          subplot(1,2,1); plot(j*ones(j,1), node_price(1:j), 'og'); set(gca,'fontsize',14);
          if j == 1
              leg = legend('Strike Price','location','best'); set(leg,'fontsize',18);
              xlabel('Number of steps','fontsize',18); ylabel('Underlting asset price','fontsize',18);
              title('Binomial tree plot of underlying asset price','fontsize',20);
          end
          hold on;
          subplot(1,2,2); plot(j*ones(j,1), put_price(1:j), '+b'); set(gca,'fontsize',14);
          xlabel('Number of steps','fontsize',18); ylabel('Option price','fontsize',18);
          title('Binomial tree plot of option price','fontsize',20); hold on;
      end
  end
  put_value = put_price(1);
end

function [call_value, exercise_max, binomial_max, exercise_time] = American_call(K,T,N,Vol,S,r)
  if T > 365
      N = floor(N*1.8);
  end
  dT = T/N/365;
  up = exp(Vol*sqrt(dT));
  down = exp(-Vol*sqrt(dT));
  up_prob = (exp(r*dT)-down) / (up-down);
  down_prob = (up-exp(r*dT)) / (up-down);
  depth = N+1;
  node_price = zeros(depth,1);
  call_price = zeros(depth,1);
  exercise_max = 0; binomial_max = 0; exercise_time = 0; whether_exercise = 0;
  for j = depth:(-1):1
      for i = 1:1:j
          node_price(i) = S*up^(j-i)*down^(i-1);
          if j == depth
              call_price(i) = max(node_price(i)-K, 0);
          else
              call_price(i) = max((node_price(i)-K), ...
                                  exp(-r*dT)*(up_prob*node_price(i)+down_prob*node_price(i+1)));
              if call_price(i) == node_price(i)-K        % analysis of early exercise
                  exercise_max = exercise_max + 1;
                  if whether_exercise == 0
                      exercise_time = j-1;
                      whether_exercise = whether_exercise+1;
                  end
              else
                  binomial_max = binomial_max + 1;
              end
          end
      end
      whether_exercise = 0;
  end
  call_value = call_price(1);
end

%% Algorithm of binomial tree model for convergence test of American option pricing
function [value_col] = Convergence_test(K,T,Vol,S,r,poc)
  for num_step = 1:1:150
      dT = T/num_step/365;
      up = exp(Vol*sqrt(dT));
      down = exp(-Vol*sqrt(dT));
      up_prob = (exp(r*dT)-down) / (up-down);
      down_prob = (up-exp(r*dT)) / (up-down);
      depth = num_step+1;
      node_price = zeros(depth,1);
      opt_price = zeros(depth,1);
      for j = depth:(-1):1
          for i = 1:1:j
              node_price(i) = S*up^(j-i)*down^(i-1);
              if j == depth
                  if poc == 1
                      opt_price(i) = max(K-node_price(i), 0);
                  else
                      opt_price(i) = max(node_price(i)-K, 0);
                  end
              else
                  if poc == 1
                      opt_price(i) = max((K-node_price(i)), ...
                                         exp(-r*dT)*(up_prob*node_price(i)+down_prob*node_price(i+1)));
                  else
                      opt_price(i) = max((node_price(i)-K), ...
                                         exp(-r*dT)*(up_prob*node_price(i)+down_prob*node_price(i+1)));
                  end
              end
          end
      end
      value_col(num_step,1) = opt_price(1);
  end
end