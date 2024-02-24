clc;
clear all;
close all;



function independentSet()
%independent set
nodes = randi(50,1,10)
%nodes = [ 42    41     4    20    27    21    33    32    15    22]
n = length(nodes);
c = zeros(n,1);
s = cell(n,1);
c(1) = nodes(1);
s{1} = 'current';
for i = 2:n
    q1 = c(i-1);
    if (i == 2)
        q2 = nodes(i);
    else
        q2 = nodes(i) + c(i - 2);
    end
    if (q1 > q2)
        s{i} = 'back';
    else
        s{i} = 'current';
    end
    c(i) = max(q1,q2);
end
s = printSet(nodes,s,n)
s = strsplit(s{1},{' '});
soma = 0;
for i = 2:length(s)
    id = str2num(s{i});
    soma = soma + nodes(id);
end
soma
c(end)
end
function s = printSet(nodes,backtracking,i)
    if (i < 1)
        s = [];
        return
    end
    if strcmp(backtracking{i},'back')
        s = strcat(printSet(nodes,backtracking,i-1),{' '});
    else
        s = strcat(printSet(nodes,backtracking,i-2),{' '},int2str(i));
    end
end

function LCS()
X = {'A','B','C','B','D','A','B'};
Y = {'B','D','C','A','B','A'};

m = length(X);
n = length(Y);
c = zeros(m,n);
left = zeros(m,n);
up = zeros(m,n);
for i = 1:m
    for j = 1:n
        if X{i} == Y{j}
            if (i == 1 || j == 1)
                c(i,j) = 1;
                left(i,j) = i;
                up(i,j) = j;
            else
                c(i,j) = c(i-1,j-1) + 1;
                left(i,j) = i-1;  
                up(i,j) = j-1;
            end
        else
            if (i == 1)
                q1 = 0;
                left(i,j) = i; 
            else
                q1 = c(i-1,j);
                left(i,j) = i-1; 
            end
            if (j == 1)
                q2 = 0;
                up(i,j) = j; 
            else
                q2 = c(i,j-1);
                up(i,j) = j-1;
            end   
            c(i,j) = max(q1,q2);
        end
    end
end
c
end

function matrixChain()
p = [30,35,15,5,10,20,25];
n = length(p) - 1;
m = (-1)*ones(n,n);
s = m;
for i = 1:n
    m(i,i) = 0.0; % não tem custo a mesma multiplicação
end
for l = 2:n
    for i = 1: n-l + 1
        j = i + l - 1;
        m(i,j) = realmax;
        for k = i: j - 1
            c = m(i,k) + m(k+1,j) + p(i)*p(k+1)*p(j+1);
            if (c < m(i,j))
                m(i,j) = c;
                s(i,j) = k;
            end
        end
    end
end
printParenthesis(s,1,6);
end
function s = printParenthesis(s,i,j)
if (i == j)
    s = strcat('A','_',int2str(i),{' '});
    return
end
k = s(i,j);
s = strcat('(',printParenthesis(s,i,k),printParenthesis(s,k+1,j),')');
s
end

