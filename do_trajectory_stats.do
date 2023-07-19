/*
Code created Feb 20/23
Xiao Chen
chenxiaophd@gmail.com
*/

clear all
set more off

* clean data
cd /Users/ChenXiao/Documents/My_Documents/Trajectory_Project
* the last excel data 
import excel Salience_S51TTR_ITT+Opt_BilateralDLPFC_Addl.xlsx, sheet(1ProtocolStanBil=AllTx#) firstrow

*drop duplicate subjects and patients that are "opted for additional treatments"
drop if substr(EnrolledinStudy,1,5) == "Opted"


* drop patients marked as red, who were with a PHQ-9 score > 5 and < 10
drop if RID == "RID-16081"
drop if RID == "RID-00127"
drop if RID == "RID-00551"
drop if RID == "RID-16532"
drop if RID == "RID-00197"
drop if RID == "RID-01642"
drop if RID == "RID-01671"
drop if RID == "RID-01663"
drop if RID == "RID-01765"
drop if RID == "RID-00409"
drop if RID == "RID-01755"
drop if RID == "RID-46997"

* exclude patients with PHQ-9 scores < 10 from remaining dataset
keep if BaselinePHQ9 >= 10

* transfer PHQ-9 scores from string to number
destring W1PHQ9, replace force
destring W2PHQ9, replace force
destring W3PHQ9, replace force
destring W4PHQ9, replace force
destring W5PHQ9, replace force
destring W6PHQ9, replace force
destring W7PHQ9, replace force
destring W8PHQ9, replace force
destring W9PHQ9, replace force
destring W10PHQ9, replace force

generate wk0 = 0
generate wk1 = 1
generate wk2 = 2
generate wk3 = 3
generate wk4 = 4
generate wk5 = 5
generate wk6 = 6
generate wk7 = 7
generate wk8 = 8
generate wk9 = 9
generate wk10 = 10

/*
*get the basic stats for the final sample
summarize AgeMTDay, detail
tabulate BiologicalSexandGenderifDi, missing
*/

*save cleaned data
save "Trajectory_Data_Clean.dta", replace
* save as .csv file to do analyses in R
export delimited using Trajectory_Data_Clean.csv, replace
*/

* use clean data to fit trajectory models
clear all
set more off

use Trajectory_Data_Clean.dta


* Seven week data

* step 01: try 1-5 group, all quadratic models
* one trajectory
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (2) min(0) max(27)

* two trajectories
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (2 2) min(0) max(27)

* three trajectories
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (2 2 2) min(0) max(27)

* four trajectories 
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (2 2 2 2) min(0) max(27)

* five trajectories 
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (2 2 2 2 2) min(0) max(27)

* step 02: looping across all possible order conbinations in the 4 group model
scalar BIC = -7663 //the best BIC
local count 0
forvalues order1 = 1/3 {
	forvalues order2 = 1/3 {
		forvalues order3 = 1/3 {
			forvalues order4 = 1/3 {
        quietly traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order(`order1' `order2' `order3' `order4') min(0) max(27)
	    scalar BIC_current = e(BIC_n_subjects)
		if BIC_current > BIC {
			local count = `count' + 1
			matrix BIC_optimized`count' = BIC_current
			matrix order1_optimized`count' = `order1'
			matrix order2_optimized`count' = `order2'
			matrix order3_optimized`count' = `order3'
			matrix order4_optimized`count' = `order4'
			scalar BIC = BIC_current
				}		
			}
		}
	}
}
matrix dir

display BIC_optimized3[1,1]
display order1_optimized3[1,1]
display order2_optimized3[1,1]
display order3_optimized3[1,1]
display order4_optimized3[1,1]

* The optimized four trajectories 
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (3 3 2 3) min(0) max(27)

* plot figure
trajplot, xtitle(Week) ytitle(PHQ-9) ci


* week 10 model


* Step 01: try 1-5 groups, all quadratic models
* one trajectory
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (2) min(0) max(27)

* two trajectories
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (2 2) min(0) max(27)

* three trajectories
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (2 2 2) min(0) max(27)

* four trajectories
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (2 2 2 2) min(0) max(27)

* five trajectories
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (2 2 2 2 2) min(0) max(27)

*  Step 02: looping across all possible order conbinations
scalar BIC = -7358 //the best BIC
local count 0
forvalues order1 = 1/3 {
	forvalues order2 = 1/3 {
		forvalues order3 = 1/3 {
			forvalues order4 = 1/3 {
        quietly traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order(`order1' `order2' `order3' `order4') min(0) max(27)
	    scalar BIC_current = e(BIC_n_subjects)
		if BIC_current > BIC {
			local count = `count' + 1
			matrix BIC_optimized`count' = BIC_current
			matrix order1_optimized`count' = `order1'
			matrix order2_optimized`count' = `order2'
			matrix order3_optimized`count' = `order3'
			matrix order4_optimized`count' = `order4' 
			scalar BIC = BIC_current
				}		
			}
		}
	}
}
matrix dir

display BIC_optimized10[1,1]
display order1_optimized10[1,1]
display order2_optimized10[1,1]
display order3_optimized10[1,1]
display order4_optimized10[1,1]

* select the most optimized four group model that is based on 10 week data
traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (3 3 2 3) min(0) max(27)

* plot figure
trajplot, xtitle(Week) ytitle(PHQ-9) ci

* export all variables with group info (week 10) for R
export delimited using Trajectory_Data_Clean_withGroupInfo.csv, replace

*export group variable
keep _traj_Group
export delimited using wk10_model_group.csv, replace

* also fit the the correspondeing model using 7 week data
* use clean data to fit trajectory models
clear all
set more off

use Trajectory_Data_Clean.dta

traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (3 3 2 3) min(0) max(27)

* export group variable
keep _traj_Group
export delimited using wk7_model_group.csv, replace

* get the mean posterior probability
* week 10 
clear all
set more off

use Trajectory_Data_Clean.dta

traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9 W8PHQ9 W9PHQ9 W10PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7 wk8 wk9 wk10) order (3 3 2 3) min(0) max(27)

collapse (mean) _traj_ProbG1 _traj_ProbG2 _traj_ProbG3 _traj_ProbG4, by(_traj_Group)

* week 7
clear all
set more off

use Trajectory_Data_Clean.dta

traj, model(cnorm) var(BaselinePHQ9 W1PHQ9 W2PHQ9 W3PHQ9 W4PHQ9 W5PHQ9 W6PHQ9 W7PHQ9) indep(wk0 wk1 wk2 wk3 wk4 wk5 wk6 wk7) order (3 3 1 1) min(0) max(27)

collapse (mean) _traj_ProbG1 _traj_ProbG2 _traj_ProbG3 _traj_ProbG4, by(_traj_Group)
