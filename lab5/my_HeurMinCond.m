function [result] = my_HeurMinCond(A_inf, A_sup, iter_number)
  
%   определяем размеры данной матрицы 
m = size(A_inf,1); 
n = size(A_inf,2);   
 
%   задаём количество случайных бросаний в реализуемом алгоритме 
if(nargin >= 3)
    NN = iter_number; 
else
    NN = 10;
end

%   инициализируем угловые матрицы для A 
Matr1 = ones(m,n); 
Matr2 = ones(m,n);   
  
%   инициализируем MinCond - минимум чисел обусловленности точечных 
%   матриц, содержащихся в заданной интервальной матрице A   
MinCond = Inf; 
  
  
for k = 1:NN 

    %   случайно порождаем целочисленную матрицу EPM из нулей и 
    %   единиц, тех же размеров,  что  и  A  (интервал случайных 
    %   целочисленных значений указывается первым аргументом randi) 
    EPM = randi([0,1],m,n); 
      
    %   порождаем угловые матрицы, диагонально противоположные 
    %   друг другу, в соответствии с эвристикой "диагональных" 
    %   методов оптимизации     
    for i = 1:m
        for j = 1:n
            if EPM(i,j) == 0 
                Matr1(i,j) = A_inf(i,j); 
                Matr2(i,j) = A_sup(i,j); 
            else 
                Matr1(i,j) = A_sup(i,j);
                Matr2(i,j) = A_inf(i,j);                 
            end 
        end
    end 
    
    %   находим числа обусловленности полученных 
    %   угловых матриц, корректируем оценку минимума  
    c1 = cond(Matr1,2); 
    c2 = cond(Matr2,2); 
    if MinCond > c1 
        MinCond = c1; 
    end
    if MinCond > c2 
        MinCond = c2; 
    end
    
end

result = MinCond;

%   выводим найденный минимум чисел обусловленности 
%disp(MinCond); 
end
  