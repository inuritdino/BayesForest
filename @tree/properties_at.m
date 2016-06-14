function [x, s] = properties_at(obj,h)

[sp,r] = obj.sorted_trunk;

if not(isempty(sp))
    H = sp(:,3);
    h = h(:);

    N = size(h,1);

    x = zeros(N,2);
    s = zeros(N,1);

    for i = 1:N

        k = find(H == h(i),1);
        if isempty(k)
            j = find(H < h(i),1,'last');

            % Height h is higher than any in the stem.
            if isempty(j)

            % Height h is lower than any in the stem.
            elseif j == 1

            % Somewhere in the middle.
            else

                ratio = (h(i)-H(j))/(H(j+1)-H(j));
                s(i) = r(j)+abs(r(j)-r(j+1))*ratio;
                x(i,:) = sp(j+1,1:2)*ratio + sp(j,1:2)*(1-ratio);
            end
        else
            s(i) = r(k);
            x(i,:) = sp(k,1:2);
        end
    end
    
    x = [x, h];
else
     error('No stem in tree. Unable to proceed.')
end