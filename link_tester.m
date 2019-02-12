%% Blank Workspace

clear
close
clc

%% Initialize Options

link_length = 17.78:2.54:30.48;

ground_pos_x = 10:5:55;
ground_pos_y = 10:5:40;

angle = 0:10:350;

%% Create Cost Array

configs = {ones(2,1), ones(2,1), ones(2,1), ones(2,1), 0, strings, ones(2,1), ones(2,1), ones(2,1), ones(2,1);
    ones(2,1), ones(2,1), ones(2,1), ones(2,1), 0, strings, ones(2,1), ones(2,1), ones(2,1), ones(2,1);
    ones(2,1), ones(2,1), ones(2,1), ones(2,1), 0, strings, ones(2,1), ones(2,1), ones(2,1), ones(2,1);
    ones(2,1), ones(2,1), ones(2,1), ones(2,1), 0, strings, ones(2,1), ones(2,1), ones(2,1), ones(2,1);
    ones(2,1), ones(2,1), ones(2,1), ones(2,1), 0, strings, ones(2,1), ones(2,1), ones(2,1), ones(2,1)};

%% Loop Through Possibilites

for g = 1:length(link_length)
    for h = 1:length(link_length)
        for k = 1:length(ground_pos_x)
            for m = 1:length(ground_pos_y)
                for p = 1:length(ground_pos_x)
                    for q = 1:length(ground_pos_y)
                        
                        break_flag = 0;
                        
                        for r = 1:length(angle)
                            for s = 1:length(angle)
                        
                                pos_left = [ground_pos_x(k), ground_pos_y(m)];
                                pos_right = [ground_pos_x(p), ground_pos_y(q)];

                                length_left = link_length(g);
                                length_right = link_length(h);

                                angle_left = angle(r);
                                angle_right = angle(s);

                                if pdist2(pos_left, pos_right, 'euclidean') <= 10 ...
                                        || pos_left(1) > pos_right(1)
                                    
                                    C = 0;
                                    break_flag = 1;
                                    break;
                                    
                                else
                                    
                                    W = length_left .* (cosd(angle_left) + sind(angle_left) * 1i);
                                    W_ground = [1 0 pos_left(1); 0 1 pos_left(2); 0 0 1] * [real(W); imag(W); 1];
                                    
                                    Z = [1 0 pos_right(1) - W_ground(1); 0 1 pos_right(2) - W_ground(2); 0 0 1] * [real(length_right .* (cosd(angle_right) + sind(angle_right) * 1i)); imag(length_right .* (cosd(angle_right) + sind(angle_right) * 1i)); 1];
                                    Z_ground = [1 0 pos_right(1); 0 1 pos_right(2); 0 0 1] * [real(length_right .* (cosd(angle_right) + sind(angle_right) * 1i)); imag(length_right .* (cosd(angle_right) + sind(angle_right) * 1i)); 1];
                                    
                                    if W_ground(1) <= 0 || W_ground(1) >= ground_pos_x(end) || W_ground(2) <=0 || W_ground(2) >= ground_pos_y(end) ...
                                            || Z_ground(1) <= 0 || Z_ground(1) >= ground_pos_x(end) || Z_ground(2) <= 0 || Z_ground(2) >= ground_pos_y(end)
                                        
                                        C = 0;
                                        
                                    else                                   
                                        
                                        psy = atand(Z(2)./Z(1));
                                        
                                        if Z(2) > 0 && Z(1) < 0
                                            psy = psy + 180;
                                        elseif Z(2) < 0 && Z(1) < 0
                                            psy = psy - 180;
                                        end
                                        
                                        psy_prime = psy + 30;
                                        
                                        Z_prime = norm(Z) .* (cosd(psy_prime) + sind(psy_prime) * 1i);
                                        
                                        free_link = [1 0 pos_left(1); 0 1 pos_left(2); 0 0 1] * [real(Z_prime); imag(Z_prime); 1];
                                        
                                        D = sqrt((pos_right(1) - free_link(1)).^2 + (pos_right(2) - free_link(2)).^2);
                                        del = (1./4) .* sqrt((D + length_left + length_right) .* (D + length_left - length_right) .* (D - length_left + length_right) .* (-D + length_left + length_right));
                                        
                                        if imag(del) > 0
                                            
                                            C = 0;
                                            break_flag = 1;
                                            break;
                                            
                                        else
                                        
                                            Z_ground_pos(1) = ((free_link(1) + pos_right(1)) ./ 2) + (((pos_right(1) - free_link(1)) .* (length_left.^2 - length_right.^2)) ./ (2 .* D.^2)) + (2 .* ((free_link(2) - pos_right(2)) ./ (D.^2)) .* del);
                                            Z_ground_neg(1) = ((free_link(1) + pos_right(1)) ./ 2) + (((pos_right(1) - free_link(1)) .* (length_left.^2 - length_right.^2)) ./ (2 .* D.^2)) - (2 .* ((free_link(2) - pos_right(2)) ./ (D.^2)) .* del);

                                            Z_ground_pos(2) = ((free_link(2) + pos_right(2)) ./ 2) + (((pos_right(2) - free_link(2)) .* (length_left.^2 - length_right.^2)) ./ (2 .* D.^2)) - (2 .* ((free_link(1) - pos_right(1)) ./ (D.^2)) .* del);
                                            Z_ground_neg(2) = ((free_link(2) + pos_right(2)) ./ 2) + (((pos_right(2) - free_link(2)) .* (length_left.^2 - length_right.^2)) ./ (2 .* D.^2)) + (2 .* ((free_link(1) - pos_right(1)) ./ (D.^2)) .* del);

                                            W_ground_pos = [1 0 Z_ground_pos(1); 0 1 Z_ground_pos(2); 0 0 1] * [real(norm(Z) .* (cosd(psy_prime + 180) + sind(psy_prime + 180) * 1i)); imag(norm(Z_prime) .* (cosd(psy_prime + 180) + sind(psy_prime + 180) * 1i)); 1];
                                            W_ground_neg = [1 0 Z_ground_neg(1); 0 1 Z_ground_neg(2); 0 0 1] * [real(norm(Z) .* (cosd(psy_prime + 180) + sind(psy_prime + 180) * 1i)); imag(norm(Z_prime) .* (cosd(psy_prime + 180) + sind(psy_prime + 180) * 1i)); 1];

                                            if (W_ground_pos(1) <= 0 || W_ground_pos(1) >= ground_pos_x(end) || W_ground_pos(2) <=0 || W_ground_pos(2) >= ground_pos_y(end) ...
                                                || Z_ground_pos(1) <= 0 || Z_ground_pos(1) >= ground_pos_x(end) || Z_ground_pos(2) <= 0 || Z_ground_pos(2) >= ground_pos_y(end)) ...
                                                && (W_ground_neg(1) <= 0 || W_ground_neg(1) >= ground_pos_x(end) || W_ground_neg(2) <=0 || W_ground_neg(2) >= ground_pos_y(end) ...
                                                || Z_ground_neg(1) <= 0 || Z_ground_neg(1) >= ground_pos_x(end) || Z_ground_neg(2) <= 0 || Z_ground_neg(2) >= ground_pos_y(end))

                                            C = 0;

                                            else

                                                W_prime = [1 0 -pos_right(1); 0 1 -pos_right(2); 0 0 1] * Z_ground;

                                                C = 10 .* abs(pos_left(1) - pos_right(1)) + 1000 .* abs(sind(atand(W_prime(2)./W_prime(1)))) + 10 .* length_right;

                                                if (W_ground_pos(1) <= 0 || W_ground_pos(1) >= ground_pos_x(end) || W_ground_pos(2) <=0 || W_ground_pos(2) >= ground_pos_y(end) ...
                                                    || Z_ground_pos(1) <= 0 || Z_ground_pos(1) >= ground_pos_x(end) || Z_ground_pos(2) <= 0 || Z_ground_pos(2) >= ground_pos_y(end))

                                                saved_orient = "neg";
                                                
                                                C = C - (10 .* abs(W_ground(1) - W_ground_neg(1))) - (10 .* abs(Z_ground(1) - Z_ground_neg(1)));

                                                elseif (W_ground_neg(1) <= 0 || W_ground_neg(1) >= ground_pos_x(end) || W_ground_neg(2) <=0 || W_ground_neg(2) >= ground_pos_y(end) ...
                                                    || Z_ground_neg(1) <= 0 || Z_ground_neg(1) >= ground_pos_x(end) || Z_ground_neg(2) <= 0 || Z_ground_neg(2) >= ground_pos_y(end))

                                                saved_orient = "pos";
                                                
                                                C = C - (10 .* abs(W_ground(1) - W_ground_pos(1))) - (10 .* abs(Z_ground(1) - Z_ground_pos(1)));

                                                else

                                                    saved_orient = "both";
                                                    
                                                    C = 0;

                                                end
                                            end
                                        end
                                    end
                                end
                                
                                if C > configs{end,3}
                                    
                                    configs(end,:) = {pos_left', W_ground(1:2), Z_ground(1:2), pos_right', C, saved_orient, W_ground_pos(1:2), W_ground_neg(1:2), Z_ground_pos', Z_ground_neg'};
                                    configs = sortrows(configs, 5, 'descend');
                                end
                                    
                                
                            end
                            
                            if break_flag == 1
                                break;
                            end
                            
                        end 
                    end
                end
            end
        end
    end
end
