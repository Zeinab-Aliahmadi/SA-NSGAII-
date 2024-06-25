# SA-NSGAII
Self-adaptive Non-Dominated Sorting Genetic Algorithm II for Personalized Tip Design Planning
=======================================================================================
The provided coding is related to Self-adaptive NSGA-II base on real data of Montreal ity, Canada.
Copyright <2024> reserved by Seyed Zeinab Aliahmadi
*It is forbidden to share the content of this document either in whole or in parts without citing the following paper. The reader who wishes to print or save this document on any media must first get the permission of the author.

Cite the following paper:
Seyed Zeinab Aliahmadi, Armin Jabbarzadeh, Lucas A. Hof (2024), A Multi-objective Optimization Approach for Sustainable and Personalized Trip Planning:  A Self-adaptive Evolutionary Algorithm with Case Study, Expert Systems With Applications, First Revision
============================================================================================================================================================
Project description:
In contemporary tourism, integrating sustainability with personalized travel planning signifies a pivotal transformation. 
This paper introduces an innovative multi-objective optimization method for developing travel itineraries that are both sustainable and customized.
The model aims to offer personalized recommendations to tourists, focusing on reducing overall costs and environmental impacts while enhancing individual preferences. 
These recommendations include the selection of accommodations and rest spots, the daily order of attractions, detailed visit schedules, and favored transportation modes. 
The model considers the complexities of urban multi-modal transportation systems, including elements like traffic signals, weather conditions,
 and unique characteristics of each transport mode, ensuring tourists can visit attractions within their available time. 
To manage uncertainties in travel time, costs, and visit durations, a practical fuzzy optimization technique is utilized. 
A self-adaptive evolutionary algorithm based on the Non-Dominated Sorting Genetic Algorithm II (SA-NSGAII) is developed for efficiently solving the multi-objective optimization problem. 
The empirical case study using real data from Montreal demonstrates the model's application in sustainable and personalized trip planning. 
The results analysis shows how this approach can help tourists achieve a balance between sustainability goals and well-informed, preference-based decisions.

================================================================================================================================================================
To run the main file named 'SaNSGA2', the attached excel file named 'CaseStudy'  needs to be read at first.


Nh:   				%The total number of available hotels (H)%
Ns:     			%The total number of available scenic spots(S)%
Nr: 			        %The total number of available restaurants(R)%
Nsr:Ns+Nr                       %The total number of available scenic spots+restaurants(R+S)%
N:Nh+Ns+Nr                      %The total number of nodes(n)%
Nv:				%The total number of days included in a tour(V) %
t: 				%The total number of days included in a tour(T) %

xr: 				%The longitude of hotels, tourism attractions and restaurants%   
yr: 				%The latitude of hotels, tourism attractions and restaurants%   

TL1:				%The longitude of traffic lights%   Reference link: https://open.canada.ca/data/en/dataset/02ebdab9-cbf3-4f56-8c29-79fa0ed0ed2e%
TL2:				%The latitude of traffic lights%    Reference link: https://open.canada.ca/data/en/dataset/02ebdab9-cbf3-4f56-8c29-79fa0ed0ed2e%

NL(i,j):  			%The lnumber of traffic lights between two nodes i and j%

LTW:				%The Lower bound of time windows to visit a place (R+S members)%  
UTW:				%The Uper bound of time windows to visit a  place (R+S members)%

(OT1,OT2,OT3): 			%Three fuzzy numbers of operational time to visite a place (R+S members)%
UT:				%Global (utility) rate of each restaurant (R members)%

(RFC1,RFC2,RFC3): 		%Three fuzzy numbers of restaurant's food & service cost (R members)%

UTH:				%Global (utility) rate of each hotel (H members)%

(RCH1,RCH2,RCH3): 		%The booking cost of each available hotel  per day per person (H members)%

SCs:				%The ticket price of each listed SS  per day per person (S members)%

VC:				%The variable cost of each vehicle for 1 unit distance%

FC:				%The variable cost (ticket price) of each vehicle for 1 unit distance%


(v_star1,v_star2,v_star3): 		%Average velocity (AV) of each vehicle (km/h)%

C:					%carbon emission of each vehicle for 1 unit distance per passenger

B:					%Budget%

DT:					%Departure time of visitor from the hotel%

wt:					%Average waiting time at a traffic light%

MT:					%Maximum number of trips per day=maximum number of visited SS per day%

***************************
All data are based on real places in Montreal, Candada. The data is collected in 2024.

Seyed Zeinab Aliahmadi
PhD student at École de technologie supérieure - ÉTS Montréal
Systems Engineering Department

25th of June 2024
***************************




