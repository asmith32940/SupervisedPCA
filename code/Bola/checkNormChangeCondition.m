function [changeVal,changeFlag] =checkNormChangeCondition(wOld, wNew)

global NORM_MIN;

changeVal =norm(wOld-wNew,'fro');

changeFlag =0;
if(changeVal < NORM_MIN)
    changeFlag =1;
end
