%% 计算最大公因数函数
% matlab自带的gcd函数只能计算整数的，这里是把小数乘10的倍数乘到整数再计算
function result = GCD(array)
    n = length(array);
    tmp = array;
    exp = 0; % Record how many times the array has been multiplied by 10
    while ~is_intarray(tmp) % Whether array multiplied by 10's has become integer
        tmp = tmp * 10;
        exp = exp + 1;
    end
    result = gcd_array(uint32(tmp),n); % Calculate the GCD of the integer array which has been multiplied by 10's
    if (exp > 0) % If the array is multiplied by 10's, divide the GCD by 10's in return
        result = double(result) / 10^exp;
    end
end
     
% 判定被乘了10的数组是否已经变成整数了
function result = is_intarray(array)
    tmp = round(array);
    if (abs(array - tmp) < 10^(-10))
        result = 1;
    else
        result = 0;
    end
end

% Calculate the GCD of an integer array
function result = gcd_array(array, n)
    result = array(1);
    for i = 1:n-1
        result = gcd(result,array(i+1)); 
    end
end