// regions                                                               //
//              nam = northamerica (leader)                              //
//              adv = advanced countries: eu-15, australia & new-zeland  //
//              eas = eastern europe                                     //
//              men = mediterranean world                                //
//              lac = latin america and caribbean islands                //
//              jap = japan                                              //
//              ssa = subsaharian africa                                 //
//              rus = russian world                                      //
//              chi = chinese world                                      //
//              ind = indian world + pacific islands                     //


// population is divided into 8 generations:                             //
//      individuals aged 15 to 24 in japan(generation 1): n1t            //
//      individuals aged 25 to 34 in japan(generation 2): n2t            //
//      individuals aged 35 to 44 in japan(generation 3): n3t            //
//      individuals aged 45 to 54 in japan(generation 4): n4t            //
//      individuals aged 55 to 64 in japan(generation 5): n5t            //
//      individuals aged 65 to 74 in japan(generation 6): n6t            //
//      individuals aged 75 to 84 in japan(generation 7): n7t            //
//      individuals aged 85 to 94 in japan(generation 8): n8t            //


// deflating
//    1) variables deflated by the leader's technology (abar):           //
//            x(consumption). wag(wage),  b (pension benefits)           //
//            bs and bu (pension benefits for skilled resp.unskilled)    //
//    2) variables deflated by the the region's young generation (n1t):  //
//            lab(labor force)                                           //
//    3) variables deflated by abar (i.e.tfpnam) and by the region's n1t //
//            k(capital), gdp, z(assets), wea(wealth),
//            gov (governement's budget)                                 //


// variables                                                             //

//      xu1jap = period 1 consumption of 1 unskilled individuals in japan
//
//      xs1jap = period 1 consumption of skilled individuals in japan
//

//      xu2jap, xu3jap, xu4jap, xu5jap, xu6jap, xu7jap, xu8jap
//
//                     (period 2 to 8 consumptions for unskilled)        //
//      xs2jap, xs3jap, xs4jap, xs5jap, xs6jap, xs7jap, xs8jap
//
//                      (period 2 to 8 consumptions for skilled)         //

//      zu1jap = period 1 aggregated assets of all unskilled individuals
//
//                         in japan (period 1 assets equal 0)            //
//      zs1jap = period 1 assets of skilled individuals in japan
//
//                                  (period 1 assets equal 0)            //

//      zu2jap, zu3jap, zu4jap, zu5jap, zu6jap, zu7jap, zu8jap
//
//              (period 2 to 8 assets for unskilled)                     //
//      zs2jap, zs3jap, zs4jap, zs5jap, zs6jap, zs7jap, zs8jap
//
//              (period 2 to 8 assets for skilled)                       //
//      taujap = labor income tax
//
//      bujap = pension benefits for unskilled individuals
//
//      bsjap = pension benefits for skilled individuals
//



//      govjap = solde budgétaire
//

//      weajap = wealth in japan
//
//          (sum of assets of unskilled and skilled over one period)     //
//      wagjap = wage in japan
//
//      intratejap = (1 + interest rate) in japan
//
//      gdpjap = gdp in japan
//
//      kjap = capital in japan
//
//      labjap = labor force in japan
//
//      hjap = human capital in japan
//

// para et v.exogènes:                                                   //

//      ggnam = tfp growth rate of the leader (nam)
//
//      p2jap to p8jap = probability of being alive at time
//
//          for generation 2 (respectively generations 3,4,5,6,7,8)      //
//      phijap = proportion of skilled individuals
//
//      lams1jap = labor force partcipation rate of
//
//                    the skilled individuals of generation 1            //
//      lamu1jap = labor force partcipation rate of
//
//                    the skilled individuals of generation 1            //
//      edusjap = fraction of time spent in education for skilled
//
//      eduujap = fraction of time spent in education for unskilled
//
//      beta = preference rate of the households
//
//      alpha = technology parameter
//
//      tfpjap = total factor productivity in japan (ajap/anam)
//
//  ddyjap = public debt in japan                                        //
//  repujap = replacement rate for the unskilled                         //
//  repsjap = replacement rate for the skilled                           //
//  note:                                                                //
//              if repsjap=repujap ==> bismarckian pension system        //
//          if repsjap=repujap/humjap ==> beveridgian pension system     //




periods 50;

var

// ------------------------ advanced countries ------------------------ //

xs1adv xs2adv xs3adv xs4adv xs5adv xs6adv xs7adv xs8adv
xu1adv xu2adv xu3adv xu4adv xu5adv xu6adv xu7adv xu8adv
zs2adv zs3adv zs4adv zs5adv zs6adv zs7adv zs8adv
zu2adv zu3adv zu4adv zu5adv zu6adv zu7adv zu8adv
labadv gdpadv weaadv kadv
bsadv buadv govadv tauadv bengdpadv supadv

// ------------------------------- nam ------------------------------- //

xs1nam xs2nam xs3nam xs4nam xs5nam xs6nam xs7nam xs8nam
xu1nam xu2nam xu3nam xu4nam xu5nam xu6nam xu7nam xu8nam
zs2nam zs3nam zs4nam zs5nam zs6nam zs7nam zs8nam
zu2nam zu3nam zu4nam zu5nam zu6nam zu7nam zu8nam
labnam gdpnam  weanam knam
bsnam bunam govnam taunam bengdpnam supnam

// ------------------------------- ssa ------------------------------- //

xs1ssa xs2ssa xs3ssa xs4ssa xs5ssa xs6ssa xs7ssa xs8ssa
xu1ssa xu2ssa xu3ssa xu4ssa xu5ssa xu6ssa xu7ssa xu8ssa
zs2ssa zs3ssa zs4ssa zs5ssa zs6ssa zs7ssa zs8ssa
zu2ssa zu3ssa zu4ssa zu5ssa zu6ssa zu7ssa zu8ssa
labssa gdpssa  weassa kssa
bsssa bussa govssa taussa bengdpssa supssa

// ------------------------------- lac ------------------------------- //

xs1lac xs2lac xs3lac xs4lac xs5lac xs6lac xs7lac xs8lac
xu1lac xu2lac xu3lac xu4lac xu5lac xu6lac xu7lac xu8lac
zs2lac zs3lac zs4lac zs5lac zs6lac zs7lac zs8lac
zu2lac zu3lac zu4lac zu5lac zu6lac zu7lac zu8lac
lablac gdplac  wealac klac
bslac bulac govlac taulac bengdplac suplac

// ------------------------------- japan ------------------------------ //

xs1jap xs2jap xs3jap xs4jap xs5jap xs6jap xs7jap xs8jap
xu1jap xu2jap xu3jap xu4jap xu5jap xu6jap xu7jap xu8jap
zs2jap zs3jap zs4jap zs5jap zs6jap zs7jap zs8jap
zu2jap zu3jap zu4jap zu5jap zu6jap zu7jap zu8jap
labjap gdpjap weajap kjap
bsjap bujap govjap taujap bengdpjap supjap

// ---------------------------- rusland ------------------------------- //

xs1rus xs2rus xs3rus xs4rus xs5rus xs6rus xs7rus xs8rus
xu1rus xu2rus xu3rus xu4rus xu5rus xu6rus xu7rus xu8rus
zs2rus zs3rus zs4rus zs5rus zs6rus zs7rus zs8rus
zu2rus zu3rus zu4rus zu5rus zu6rus zu7rus zu8rus
labrus gdprus wearus krus
bsrus burus govrus taurus bengdprus suprus

// ------------------------------ mena -------------------------------- //

xs1men xs2men xs3men xs4men xs5men xs6men xs7men xs8men
xu1men xu2men xu3men xu4men xu5men xu6men xu7men xu8men
zs2men zs3men zs4men zs5men zs6men zs7men zs8men
zu2men zu3men zu4men zu5men zu6men zu7men zu8men
labmen gdpmen  weamen kmen
bsmen bumen govmen taumen bengdpmen supmen

// -------------------------- eastern europe -------------------------- //

xs1eas xs2eas xs3eas xs4eas xs5eas xs6eas xs7eas xs8eas
xu1eas xu2eas xu3eas xu4eas xu5eas xu6eas xu7eas xu8eas
zs2eas zs3eas zs4eas zs5eas zs6eas zs7eas zs8eas
zu2eas zu3eas zu4eas zu5eas zu6eas zu7eas zu8eas
labeas gdpeas weaeas keas
bseas bueas goveas taueas bengdpeas supeas

// ------------------------------ china ------------------------------- //

xs1chi xs2chi xs3chi xs4chi xs5chi xs6chi xs7chi xs8chi
xu1chi xu2chi xu3chi xu4chi xu5chi xu6chi xu7chi xu8chi
zs2chi zs3chi zs4chi zs5chi zs6chi zs7chi zs8chi
zu2chi zu3chi zu4chi zu5chi zu6chi zu7chi zu8chi
labchi gdpchi weachi kchi
bschi buchi govchi tauchi bengdpchi supchi

// ------------------------------ india ------------------------------- //

xs1ind xs2ind xs3ind xs4ind xs5ind xs6ind xs7ind xs8ind
xu1ind xu2ind xu3ind xu4ind xu5ind xu6ind xu7ind xu8ind
zs2ind zs3ind zs4ind zs5ind zs6ind zs7ind zs8ind
zu2ind zu3ind zu4ind zu5ind zu6ind zu7ind zu8ind
labind gdpind weaind kind
bsind buind govind tauind bengdpind supind

trgdpadv trgdpnam trgdpssa trgdplac trgdpjap trgdprus trgdpmen trgdpeas
trgdpchi
trgdpind

lablabsshadv lablabsshnam lablabsshjap lablabsshssa lablabsshlac
lablabssheas lablabsshmen lablabsshrus lablabsshchi lablabsshind

conspubgdpadv conspubgdpnam conspubgdpjap conspubgdpeas conspubgdpmen
conspubgdprus conspubgdplac conspubgdpssa conspubgdpchi conspubgdpind

labsadv labsnam labsjap labsssa labslac labseas labsmen labsrus labschi
labsind

labuadv labunam labujap labussa labulac labueas labumen laburus labuchi
labuind

labdduadv labddunam labddujap labddussa labddulac labddueas
labddumen labddurus labdduchi labdduind

wagsadv wagsnam wagsssa wagslac wagsjap wagseas wagsmen wagsrus wagschi
wagsind

waguadv wagunam wagussa wagulac wagujap wagueas wagumen wagurus waguchi
waguind

wagddsadv wagddsnam wagddsssa wagddslac wagddsjap
wagddseas wagddsmen wagddsrus wagddschi wagddsind

wagdduadv wagddunam wagddussa wagddulac wagddujap
wagddueas wagddumen wagddurus wagdduchi wagdduind

thetaadv thetanam thetassa thetalac thetajap thetamen thetaeas thetarus
thetachi
thetaind

pissa pilac pimen pieas pirus pichi piind

taxesgdpadv taxesgdpchi taxesgdpeas taxesgdpind taxesgdpjap
taxesgdplac taxesgdpmen taxesgdpnam taxesgdprus taxesgdpssa

growthadv growthnam growthjap growthssa growthlac
growthmen growthrus growtheas growthchi growthind

ownadv ownnam ownjap ownssa ownlac ownmen ownrus owneas ownchi ownind
faadv fanam fajap falac faeas famen farus fachi faind fassa
caadv canam cajap calac caeas camen carus cachi caind cassa

fasum casum ownsum

zzdebtadv zzdebtnam zzdebtjap zzdebtlac zzdebteas zzdebtmen zzdebtrus
zzdebtchi zzdebtind
zzdebtssa

zzown2adv zzown2nam zzown2jap zzown2ssa zzown2lac zzown2men zzown2rus
zzown2eas zzown2chi
zzown2ind
zzfa2adv zzfa2nam zzfa2jap zzfa2lac zzfa2eas zzfa2men zzfa2rus zzfa2chi
zzfa2ind
zzfa2ssa
zzca2adv zzca2nam zzca2jap zzca2lac zzca2eas zzca2men zzca2rus zzca2chi
zzca2ind
zzca2ssa

zzown2sum zzfa2sum zzca2sum zzdebtsum

zzvardebtadv zzvardebtnam zzvardebtjap zzvardebtlac zzvardebteas zzvardebtmen
zzvardebtrus zzvardebtchi zzvardebtind
zzvardebtssa
zzvarkadv zzvarknam zzvarkjap zzvarklac zzvarkeas zzvarkmen zzvarkrus
zzvarkchi zzvarkind
zzvarkssa
zzvarweaadv zzvarweanam zzvarweajap zzvarwealac zzvarweaeas zzvarweamen
zzvarwearus zzvarweachi zzvarweaind
zzvarweassa

zzznewcaadv zzznewcanam zzznewcajap zzznewcalac zzznewcaeas zzznewcamen
zzznewcarus zzznewcachi zzznewcaind
zzznewcassa
zzznewcasum

intrate

ratiogdpadv
ratiogdpssa
ratiogdplac
ratiogdpjap
ratiogdprus
ratiogdpmen
ratiogdpeas
ratiogdpchi
ratiogdpind

urnam
uradv
ureas
urrus
urmen
urssa
urlac
urjap
urind
urchi

hadv hnam hssa hlac hjap hrus hmen heas hchi hind

;

varexo

// ------------------------ advanced countries ------------------------ //

p2adv p3adv p4adv p5adv p6adv p7adv p8adv
edusadv eduuadv
lams1adv lams2adv lams3adv lams4adv lams5adv lams6adv lams7adv lams8adv
lamu1adv lamu2adv lamu3adv lamu4adv lamu5adv lamu6adv lamu7adv lamu8adv
ddyadv conspubadv

// ------------------------------- nam ------------------------------- //

p2nam p3nam p4nam p5nam p6nam p7nam p8nam
edusnam eduunam
lams1nam lams2nam lams3nam lams4nam lams5nam lams6nam lams7nam lams8nam
lamu1nam lamu2nam lamu3nam lamu4nam lamu5nam lamu6nam lamu7nam lamu8nam
ddynam conspubnam

// ------------------------------- ssa ------------------------------- //

p2ssa p3ssa p4ssa p5ssa p6ssa p7ssa p8ssa
edusssa eduussa
lams1ssa lams2ssa lams3ssa lams4ssa lams5ssa lams6ssa lams7ssa lams8ssa
lamu1ssa lamu2ssa lamu3ssa lamu4ssa lamu5ssa lamu6ssa lamu7ssa lamu8ssa
ddyssa conspubssa

// ------------------------------- lac ------------------------------- //

p2lac p3lac p4lac p5lac p6lac p7lac p8lac
eduslac eduulac
lams1lac lams2lac lams3lac lams4lac lams5lac lams6lac lams7lac lams8lac
lamu1lac lamu2lac lamu3lac lamu4lac lamu5lac lamu6lac lamu7lac lamu8lac
ddylac conspublac

// ------------------------------ japan ------------------------------- //

p2jap p3jap p4jap p5jap p6jap p7jap p8jap
edusjap eduujap
lams1jap lams2jap lams3jap lams4jap lams5jap lams6jap lams7jap lams8jap
lamu1jap lamu2jap lamu3jap lamu4jap lamu5jap lamu6jap lamu7jap lamu8jap
ddyjap conspubjap

// ---------------------------- rusland ------------------------------- //

p2rus p3rus p4rus p5rus p6rus p7rus p8rus
edusrus eduurus
lams1rus lams2rus lams3rus lams4rus lams5rus lams6rus lams7rus lams8rus
lamu1rus lamu2rus lamu3rus lamu4rus lamu5rus lamu6rus lamu7rus lamu8rus
ddyrus conspubrus

// ------------------------------ mena -------------------------------- //

p2men p3men p4men p5men p6men p7men p8men
edusmen eduumen
lams1men lams2men lams3men lams4men lams5men lams6men lams7men lams8men
lamu1men lamu2men lamu3men lamu4men lamu5men lamu6men lamu7men lamu8men
ddymen conspubmen

// ------------------------ eastern countries ------------------------- //

p2eas p3eas p4eas p5eas p6eas p7eas p8eas
eduseas eduueas
lams1eas lams2eas lams3eas lams4eas lams5eas lams6eas lams7eas lams8eas
lamu1eas lamu2eas lamu3eas lamu4eas lamu5eas lamu6eas lamu7eas lamu8eas
ddyeas conspubeas

// ------------------------------ china ------------------------------- //

p2chi p3chi p4chi p5chi p6chi p7chi p8chi
eduschi eduuchi
lams1chi lams2chi lams3chi lams4chi lams5chi lams6chi lams7chi lams8chi
lamu1chi lamu2chi lamu3chi lamu4chi lamu5chi lamu6chi lamu7chi lamu8chi
ddychi conspubchi

// ------------------------------ india ------------------------------- //

p2ind p3ind p4ind p5ind p6ind p7ind p8ind
edusind eduuind
lams1ind lams2ind lams3ind lams4ind lams5ind lams6ind lams7ind lams8ind
lamu1ind lamu2ind lamu3ind lamu4ind lamu5ind lamu6ind lamu7ind lamu8ind
ddyind conspubind

ggnam

repnam repadv repjap repssa replac repmen repeas reprus repchi repind

mmnam mmadv mmeas mmrus mmmen mmssa mmlac mmjap mmind mmchi

nnadvnam nnssanam nnlacnam nnjapnam nnrusnam nnmennam nneasnam nnchinam
nnindnam

phinam phiadv phijap phieas phirus phimen phissa philac phichi phiind

tcadv tcnam tcjap tclac tcssa tcmen tceas tcrus tcchi tcind

psinam psiadv psijap psieas psimen psirus psilac psissa psichi psiind

rankssa ranklac rankmen rankeas rankrus rankchi rankind

tfpadv
tfpssa
tfplac
tfpjap
tfprus
tfpmen
tfpeas
tfpchi
tfpind

etaadv etanam etajap etaeas etamen etarus etalac etassa etachi etaind

aaadv
aanam
aassa
aalac
aajap
aaeas
aamen
aarus
aachi
aaind

pimax
;

parameters
beta alpha delta
rhonam rhoadv rhojap rhossa rholac rhomen rhoeas rhorus rhochi rhoind
tr1unam tr2unam tr3unam tr4unam tr5unam tr6unam tr7unam tr8unam
tr1snam tr2snam tr3snam tr4snam tr5snam tr6snam tr7snam tr8snam
tr1uadv tr2uadv tr3uadv tr4uadv tr5uadv tr6uadv tr7uadv tr8uadv
tr1sadv tr2sadv tr3sadv tr4sadv tr5sadv tr6sadv tr7sadv tr8sadv
tr1ujap tr2ujap tr3ujap tr4ujap tr5ujap tr6ujap tr7ujap tr8ujap
tr1sjap tr2sjap tr3sjap tr4sjap tr5sjap tr6sjap tr7sjap tr8sjap
tr1ueas tr2ueas tr3ueas tr4ueas tr5ueas tr6ueas tr7ueas tr8ueas
tr1seas tr2seas tr3seas tr4seas tr5seas tr6seas tr7seas tr8seas
tr1ussa tr2ussa tr3ussa tr4ussa tr5ussa tr6ussa tr7ussa tr8ussa
tr1sssa tr2sssa tr3sssa tr4sssa tr5sssa tr6sssa tr7sssa tr8sssa
tr1ulac tr2ulac tr3ulac tr4ulac tr5ulac tr6ulac tr7ulac tr8ulac
tr1slac tr2slac tr3slac tr4slac tr5slac tr6slac tr7slac tr8slac
tr1umen tr2umen tr3umen tr4umen tr5umen tr6umen tr7umen tr8umen
tr1smen tr2smen tr3smen tr4smen tr5smen tr6smen tr7smen tr8smen
tr1urus tr2urus tr3urus tr4urus tr5urus tr6urus tr7urus tr8urus
tr1srus tr2srus tr3srus tr4srus tr5srus tr6srus tr7srus tr8srus
tr1uchi tr2uchi tr3uchi tr4uchi tr5uchi tr6uchi tr7uchi tr8uchi
tr1schi tr2schi tr3schi tr4schi tr5schi tr6schi tr7schi tr8schi
tr1uind tr2uind tr3uind tr4uind tr5uind tr6uind tr7uind tr8uind
tr1sind tr2sind tr3sind tr4sind tr5sind tr6sind tr7sind tr8sind

sigma
;

sigma=0.285714286;
alpha=0.33;
beta=1;
delta=0.4;
rhonam=0.2;
rhoadv=0.6;
rhojap=0.8;
rhossa=0;
rholac=0;
rhomen=0;
rhoeas=0;
rhorus=0;
rhochi=0;
rhoind=0;

tr1unam=        0.07    ;
tr2unam=        0.12    ;
tr3unam=        0.12    ;
tr4unam=        0.12    ;
tr5unam=        0.12    ;
tr6unam=        0.2     ;
tr7unam=        0.25    ;
tr8unam=        0.25    ;
tr1snam=        0.005   ;
tr2snam=        0.01    ;
tr3snam=        0.01    ;
tr4snam=        0.01    ;
tr5snam=        0.01    ;
tr6snam=        0.06    ;
tr7snam=        0.12    ;
tr8snam=        0.12    ;

tr1uadv=        0.1     ;
tr2uadv=        0.18    ;
tr3uadv=        0.18    ;
tr4uadv=        0.18    ;
tr5uadv=        0.18    ;
tr6uadv=        0.2     ;
tr7uadv=        0.22    ;
tr8uadv=        0.22    ;
tr1sadv=        0.03    ;
tr2sadv=        0.06    ;
tr3sadv=        0.06    ;
tr4sadv=        0.06    ;
tr5sadv=        0.06    ;
tr6sadv=        0.1     ;
tr7sadv=        0.15    ;
tr8sadv=        0.15    ;

tr1ueas=        0.007   ;
tr2ueas=        0.012   ;
tr3ueas=        0.012   ;
tr4ueas=        0.012   ;
tr5ueas=        0.012   ;
tr6ueas=        0.02    ;
tr7ueas=        0.025   ;
tr8ueas=        0.025   ;
tr1seas=        0.0005  ;
tr2seas=        0.001   ;
tr3seas=        0.001   ;
tr4seas=        0.001   ;
tr5seas=        0.001   ;
tr6seas=        0.006   ;
tr7seas=        0.012   ;
tr8seas=        0.012   ;

tr1ulac=        0.007   ;
tr2ulac=        0.012   ;
tr3ulac=        0.012   ;
tr4ulac=        0.012   ;
tr5ulac=        0.012   ;
tr6ulac=        0.02    ;
tr7ulac=        0.025   ;
tr8ulac=        0.025   ;
tr1slac=        0.0005  ;
tr2slac=        0.001   ;
tr3slac=        0.001   ;
tr4slac=        0.001   ;
tr5slac=        0.001   ;
tr6slac=        0.006   ;
tr7slac=        0.012   ;
tr8slac=        0.012   ;

tr1ussa=        0.007   ;
tr2ussa=        0.012   ;
tr3ussa=        0.012   ;
tr4ussa=        0.012   ;
tr5ussa=        0.012   ;
tr6ussa=        0.02    ;
tr7ussa=        0.025   ;
tr8ussa=        0.025   ;
tr1sssa=        0.0005  ;
tr2sssa=        0.001   ;
tr3sssa=        0.001   ;
tr4sssa=        0.001   ;
tr5sssa=        0.001   ;
tr6sssa=        0.006   ;
tr7sssa=        0.012   ;
tr8sssa=        0.012   ;

tr1uchi=        0.007   ;
tr2uchi=        0.012   ;
tr3uchi=        0.012   ;
tr4uchi=        0.012   ;
tr5uchi=        0.012   ;
tr6uchi=        0.02    ;
tr7uchi=        0.025   ;
tr8uchi=        0.025   ;
tr1schi=        0.0005  ;
tr2schi=        0.001   ;
tr3schi=        0.001   ;
tr4schi=        0.001   ;
tr5schi=        0.001   ;
tr6schi=        0.006   ;
tr7schi=        0.012   ;
tr8schi=        0.012   ;

tr1uind=        0.007   ;
tr2uind=        0.012   ;
tr3uind=        0.012   ;
tr4uind=        0.012   ;
tr5uind=        0.012   ;
tr6uind=        0.02    ;
tr7uind=        0.025   ;
tr8uind=        0.025   ;
tr1sind=        0.0005  ;
tr2sind=        0.001   ;
tr3sind=        0.001   ;
tr4sind=        0.001   ;
tr5sind=        0.001   ;
tr6sind=        0.006   ;
tr7sind=        0.012   ;
tr8sind=        0.012   ;

tr1umen=        0.007   ;
tr2umen=        0.012   ;
tr3umen=        0.012   ;
tr4umen=        0.012   ;
tr5umen=        0.012   ;
tr6umen=        0.02    ;
tr7umen=        0.025   ;
tr8umen=        0.025   ;
tr1smen=        0.0005  ;
tr2smen=        0.001   ;
tr3smen=        0.001   ;
tr4smen=        0.001   ;
tr5smen=        0.001   ;
tr6smen=        0.006   ;
tr7smen=        0.012   ;
tr8smen=        0.012   ;

tr1urus=        0.007   ;
tr2urus=        0.012   ;
tr3urus=        0.012   ;
tr4urus=        0.012   ;
tr5urus=        0.012   ;
tr6urus=        0.02    ;
tr7urus=        0.025   ;
tr8urus=        0.025   ;
tr1srus=        0.0005  ;
tr2srus=        0.001   ;
tr3srus=        0.001   ;
tr4srus=        0.001   ;
tr5srus=        0.001   ;
tr6srus=        0.006   ;
tr7srus=        0.012   ;
tr8srus=        0.012   ;

tr1ujap=        0.1     ;
tr2ujap=        0.18    ;
tr3ujap=        0.18    ;
tr4ujap=        0.18    ;
tr5ujap=        0.18    ;
tr6ujap=        0.2     ;
tr7ujap=        0.22    ;
tr8ujap=        0.22    ;
tr1sjap=        0.03    ;
tr2sjap=        0.06    ;
tr3sjap=        0.06    ;
tr4sjap=        0.06    ;
tr5sjap=        0.06    ;
tr6sjap=        0.1     ;
tr7sjap=        0.15    ;
tr8sjap=        0.15    ;


model(sparse_dll,GCC_Compiler, cutoff = 1.0e-17);

// -------------------------------- open ------------------------------- //

knam+kadv*nnadvnam+kssa*nnssanam+klac*nnlacnam
+kjap*nnjapnam+krus*nnrusnam+kmen*nnmennam
+keas*nneasnam+kchi*nnchinam+kind*nnindnam
+ddynam*gdpnam
+ddyadv*gdpadv*nnadvnam+ddyssa*gdpssa*nnssanam+ddylac*gdplac*nnlacnam
+ddyjap*gdpjap*nnjapnam+ddyrus*gdprus*nnrusnam+ddymen*gdpmen*nnmennam
+ddyeas*gdpeas*nneasnam+ddychi*gdpchi*nnchinam+ddyind*gdpind*nnindnam
=weanam
+weaadv*nnadvnam+weassa*nnssanam+wealac*nnlacnam
+weajap*nnjapnam+wearus*nnrusnam+weamen*nnmennam
+weaeas*nneasnam+weachi*nnchinam+weaind*nnindnam;

intrate=1+alpha*knam^(alpha-1)*labnam^(1-alpha)-delta;
knam^(alpha-1)*labnam^(1-alpha)=kadv^(alpha-1)*(tfpadv*labadv)^(1-alpha);
knam^(alpha-1)*labnam^(1-alpha)=kjap^(alpha-1)*(tfpjap*labjap)^(1-alpha);

intrate*(1+pissa)=(1+alpha*kssa^(alpha-1)*(tfpssa*labssa)^(1-alpha)-delta);
intrate*(1+pilac)=(1+alpha*klac^(alpha-1)*(tfplac*lablac)^(1-alpha)-delta);
intrate*(1+pirus)=(1+alpha*krus^(alpha-1)*(tfprus*labrus)^(1-alpha)-delta);
intrate*(1+pimen)=(1+alpha*kmen^(alpha-1)*(tfpmen*labmen)^(1-alpha)-delta);
intrate*(1+pieas)=(1+alpha*keas^(alpha-1)*(tfpeas*labeas)^(1-alpha)-delta);
intrate*(1+pichi)=(1+alpha*kchi^(alpha-1)*(tfpchi*labchi)^(1-alpha)-delta);
intrate*(1+piind)=(1+alpha*kind^(alpha-1)*(tfpind*labind)^(1-alpha)-delta);

ddyadv(+1)*gdpadv(+1)=intrate*ddyadv*gdpadv/(ggnam(+1)*mmadv)-govadv/(ggnam(+1)*mmadv);
ddynam(+1)*gdpnam(+1)=intrate*ddynam*gdpnam/(ggnam(+1)*mmnam)-govnam/(ggnam(+1)*mmnam);
ddyssa(+1)*gdpssa(+1)=intrate*ddyssa*gdpssa/(ggnam(+1)*mmssa)-govssa/(ggnam(+1)*mmssa);
ddylac(+1)*gdplac(+1)=intrate*ddylac*gdplac/(ggnam(+1)*mmlac)-govlac/(ggnam(+1)*mmlac);
ddyjap(+1)*gdpjap(+1)=intrate*ddyjap*gdpjap/(ggnam(+1)*mmjap)-govjap/(ggnam(+1)*mmjap);
ddyrus(+1)*gdprus(+1)=intrate*ddyrus*gdprus/(ggnam(+1)*mmrus)-govrus/(ggnam(+1)*mmrus);
ddymen(+1)*gdpmen(+1)=intrate*ddymen*gdpmen/(ggnam(+1)*mmmen)-govmen/(ggnam(+1)*mmmen);
ddyeas(+1)*gdpeas(+1)=intrate*ddyeas*gdpeas/(ggnam(+1)*mmeas)-goveas/(ggnam(+1)*mmeas);
ddychi(+1)*gdpchi(+1)=intrate*ddychi*gdpchi/(ggnam(+1)*mmchi)-govchi/(ggnam(+1)*mmchi);
ddyind(+1)*gdpind(+1)=intrate*ddyind*gdpind/(ggnam(+1)*mmind)-govind/(ggnam(+1)*mmind);

bsadv=repadv*(rhoadv*wagddsadv+(1-rhoadv)*(1-uradv)*wagdduadv);
buadv=repadv*(1-uradv)*wagdduadv;
bsnam=repnam*(rhonam*wagddsnam+(1-rhonam)*(1-urnam)*wagddunam);
bunam=repnam*(1-urnam)*wagddunam;
bsssa=repssa*(rhossa*wagddsssa+(1-rhossa)*(1-urssa)*wagddussa);
bussa=repssa*(1-urssa)*wagddussa;
bslac=replac*(rholac*wagddslac+(1-rholac)*(1-urlac)*wagddulac);
bulac=replac*(1-urlac)*wagddulac;
bsjap=repjap*(rhojap*wagddsjap+(1-rhojap)*(1-urjap)*wagddujap);
bujap=repjap*(1-urjap)*wagddujap;
bsrus=reprus*(rhorus*wagddsrus+(1-rhorus)*(1-urrus)*wagddurus);
burus=reprus*(1-urrus)*wagddurus;
bsmen=repmen*(rhomen*wagddsmen+(1-rhomen)*(1-urmen)*wagddumen);
bumen=repmen*(1-urmen)*wagddumen;
bseas=repeas*(rhoeas*wagddseas+(1-rhoeas)*(1-ureas)*wagddueas);
bueas=repeas*(1-ureas)*wagddueas;
bschi=repchi*(rhochi*wagddschi+(1-rhochi)*(1-urchi)*wagdduchi);
buchi=repchi*(1-urchi)*wagdduchi;
bsind=repind*(rhoind*wagddsind+(1-rhoind)*(1-urind)*wagdduind);
buind=repind*(1-urind)*wagdduind;

// ----------------------------- pi ----------------------------------- //

pissa=rankssa/7*pimax;
pilac=ranklac/7*pimax;
pimen=rankmen/7*pimax;
pieas=rankeas/7*pimax;
pirus=rankrus/7*pimax;
pichi=rankchi/7*pimax;
piind=rankind/7*pimax;

// ----------------------------- labs and labu ------------------------ //

labsadv=(1-edusadv)*phiadv*lams1adv
+phiadv(-1)*lams2adv*p2adv/mmadv(-1)
+phiadv(-2)*lams3adv*p3adv/(mmadv(-1)*mmadv(-2))
+phiadv(-3)*lams4adv*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+phiadv(-4)*lams5adv*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+phiadv(-5)*lams6adv*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+phiadv(-6)*lams7adv*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+phiadv(-7)*lams8adv*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7));

labuadv=(1-eduuadv)*(1-phiadv)*lamu1adv
+(1-phiadv(-1))*lamu2adv*p2adv/mmadv(-1)
+(1-phiadv(-2))*lamu3adv*p3adv/(mmadv(-1)*mmadv(-2))
+(1-phiadv(-3))*lamu4adv*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+(1-phiadv(-4))*lamu5adv*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+(1-phiadv(-5))*lamu6adv*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+(1-phiadv(-6))*lamu7adv*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+(1-phiadv(-7))*lamu8adv*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7));

labsnam=(1-edusnam)*phinam*lams1nam
+phinam(-1)*lams2nam*p2nam/mmnam(-1)
+phinam(-2)*lams3nam*p3nam/(mmnam(-1)*mmnam(-2))
+phinam(-3)*lams4nam*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+phinam(-4)*lams5nam*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+phinam(-5)*lams6nam*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+phinam(-6)*lams7nam*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+phinam(-7)*lams8nam*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7));

labunam=(1-eduunam)*(1-phinam)*lamu1nam
+(1-phinam(-1))*lamu2nam*p2nam/mmnam(-1)
+(1-phinam(-2))*lamu3nam*p3nam/(mmnam(-1)*mmnam(-2))
+(1-phinam(-3))*lamu4nam*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+(1-phinam(-4))*lamu5nam*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+(1-phinam(-5))*lamu6nam*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+(1-phinam(-6))*lamu7nam*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+(1-phinam(-7))*lamu8nam*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7));

labsssa=(1-edusssa)*phissa*lams1ssa
+phissa(-1)*lams2ssa*p2ssa/mmssa(-1)
+phissa(-2)*lams3ssa*p3ssa/(mmssa(-1)*mmssa(-2))
+phissa(-3)*lams4ssa*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+phissa(-4)*lams5ssa*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+phissa(-5)*lams6ssa*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+phissa(-6)*lams7ssa*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+phissa(-7)*lams8ssa*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7));

labussa=(1-eduussa)*(1-phissa)*lamu1ssa
+(1-phissa(-1))*lamu2ssa*p2ssa/mmssa(-1)
+(1-phissa(-2))*lamu3ssa*p3ssa/(mmssa(-1)*mmssa(-2))
+(1-phissa(-3))*lamu4ssa*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+(1-phissa(-4))*lamu5ssa*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+(1-phissa(-5))*lamu6ssa*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+(1-phissa(-6))*lamu7ssa*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+(1-phissa(-7))*lamu8ssa*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7));

labslac=(1-eduslac)*philac*lams1lac
+philac(-1)*lams2lac*p2lac/mmlac(-1)
+philac(-2)*lams3lac*p3lac/(mmlac(-1)*mmlac(-2))
+philac(-3)*lams4lac*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+philac(-4)*lams5lac*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+philac(-5)*lams6lac*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+philac(-6)*lams7lac*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+philac(-7)*lams8lac*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7));

labulac=(1-eduulac)*(1-philac)*lamu1lac
+(1-philac(-1))*lamu2lac*p2lac/mmlac(-1)
+(1-philac(-2))*lamu3lac*p3lac/(mmlac(-1)*mmlac(-2))
+(1-philac(-3))*lamu4lac*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+(1-philac(-4))*lamu5lac*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+(1-philac(-5))*lamu6lac*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+(1-philac(-6))*lamu7lac*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+(1-philac(-7))*lamu8lac*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7));

labsjap=(1-edusjap)*phijap*lams1jap
+phijap(-1)*lams2jap*p2jap/mmjap(-1)
+phijap(-2)*lams3jap*p3jap/(mmjap(-1)*mmjap(-2))
+phijap(-3)*lams4jap*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+phijap(-4)*lams5jap*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+phijap(-5)*lams6jap*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+phijap(-6)*lams7jap*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+phijap(-7)*lams8jap*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7));

labujap=(1-eduujap)*(1-phijap)*lamu1jap
+(1-phijap(-1))*lamu2jap*p2jap/mmjap(-1)
+(1-phijap(-2))*lamu3jap*p3jap/(mmjap(-1)*mmjap(-2))
+(1-phijap(-3))*lamu4jap*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+(1-phijap(-4))*lamu5jap*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+(1-phijap(-5))*lamu6jap*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+(1-phijap(-6))*lamu7jap*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+(1-phijap(-7))*lamu8jap*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7));

labseas=(1-eduseas)*phieas*lams1eas
+phieas(-1)*lams2eas*p2eas/mmeas(-1)
+phieas(-2)*lams3eas*p3eas/(mmeas(-1)*mmeas(-2))
+phieas(-3)*lams4eas*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+phieas(-4)*lams5eas*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+phieas(-5)*lams6eas*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+phieas(-6)*lams7eas*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+phieas(-7)*lams8eas*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7));

labueas=(1-eduueas)*(1-phieas)*lamu1eas
+(1-phieas(-1))*lamu2eas*p2eas/mmeas(-1)
+(1-phieas(-2))*lamu3eas*p3eas/(mmeas(-1)*mmeas(-2))
+(1-phieas(-3))*lamu4eas*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+(1-phieas(-4))*lamu5eas*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+(1-phieas(-5))*lamu6eas*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+(1-phieas(-6))*lamu7eas*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+(1-phieas(-7))*lamu8eas*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7));

labsmen=(1-edusmen)*phimen*lams1men
+phimen(-1)*lams2men*p2men/mmmen(-1)
+phimen(-2)*lams3men*p3men/(mmmen(-1)*mmmen(-2))
+phimen(-3)*lams4men*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+phimen(-4)*lams5men*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+phimen(-5)*lams6men*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+phimen(-6)*lams7men*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+phimen(-7)*lams8men*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7));

labumen=(1-eduumen)*(1-phimen)*lamu1men
+(1-phimen(-1))*lamu2men*p2men/mmmen(-1)
+(1-phimen(-2))*lamu3men*p3men/(mmmen(-1)*mmmen(-2))
+(1-phimen(-3))*lamu4men*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+(1-phimen(-4))*lamu5men*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+(1-phimen(-5))*lamu6men*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+(1-phimen(-6))*lamu7men*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+(1-phimen(-7))*lamu8men*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7));

labsrus=(1-edusrus)*phirus*lams1rus
+phirus(-1)*lams2rus*p2rus/mmrus(-1)
+phirus(-2)*lams3rus*p3rus/(mmrus(-1)*mmrus(-2))
+phirus(-3)*lams4rus*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+phirus(-4)*lams5rus*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+phirus(-5)*lams6rus*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+phirus(-6)*lams7rus*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+phirus(-7)*lams8rus*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7));

laburus=(1-eduurus)*(1-phirus)*lamu1rus
+(1-phirus(-1))*lamu2rus*p2rus/mmrus(-1)
+(1-phirus(-2))*lamu3rus*p3rus/(mmrus(-1)*mmrus(-2))
+(1-phirus(-3))*lamu4rus*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+(1-phirus(-4))*lamu5rus*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+(1-phirus(-5))*lamu6rus*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+(1-phirus(-6))*lamu7rus*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+(1-phirus(-7))*lamu8rus*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7));

labschi=(1-eduschi)*phichi*lams1chi
+phichi(-1)*lams2chi*p2chi/mmchi(-1)
+phichi(-2)*lams3chi*p3chi/(mmchi(-1)*mmchi(-2))
+phichi(-3)*lams4chi*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+phichi(-4)*lams5chi*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+phichi(-5)*lams6chi*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+phichi(-6)*lams7chi*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+phichi(-7)*lams8chi*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7));

labuchi=(1-eduuchi)*(1-phichi)*lamu1chi
+(1-phichi(-1))*lamu2chi*p2chi/mmchi(-1)
+(1-phichi(-2))*lamu3chi*p3chi/(mmchi(-1)*mmchi(-2))
+(1-phichi(-3))*lamu4chi*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+(1-phichi(-4))*lamu5chi*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+(1-phichi(-5))*lamu6chi*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+(1-phichi(-6))*lamu7chi*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+(1-phichi(-7))*lamu8chi*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7));

labsind=(1-edusind)*phiind*lams1ind
+phiind(-1)*lams2ind*p2ind/mmind(-1)
+phiind(-2)*lams3ind*p3ind/(mmind(-1)*mmind(-2))
+phiind(-3)*lams4ind*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+phiind(-4)*lams5ind*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+phiind(-5)*lams6ind*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+phiind(-6)*lams7ind*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+phiind(-7)*lams8ind*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7));

labuind=(1-eduuind)*(1-phiind)*lamu1ind
+(1-phiind(-1))*lamu2ind*p2ind/mmind(-1)
+(1-phiind(-2))*lamu3ind*p3ind/(mmind(-1)*mmind(-2))
+(1-phiind(-3))*lamu4ind*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+(1-phiind(-4))*lamu5ind*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+(1-phiind(-5))*lamu6ind*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+(1-phiind(-6))*lamu7ind*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+(1-phiind(-7))*lamu8ind*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7));


// -------------------------------- aa ------------------------------- //

wagsadv/waguadv=hadv;
wagsnam/wagunam=hnam;
wagsssa/wagussa=hssa;
wagslac/wagulac=hlac;
wagsjap/wagujap=hjap;
wagseas/wagueas=heas;
wagsmen/wagumen=hmen;
wagsrus/wagurus=hrus;
wagschi/waguchi=hchi;
wagsind/waguind=hind;

wagddsadv/wagdduadv=thetaadv;
wagddsnam/wagddunam=thetanam;
wagddsssa/wagddussa=thetassa;
wagddslac/wagddulac=thetalac;
wagddsjap/wagddujap=thetajap;
wagddseas/wagddueas=thetaeas;
wagddsmen/wagddumen=thetamen;
wagddsrus/wagddurus=thetarus;
wagddschi/wagdduchi=thetachi;
wagddsind/wagdduind=thetaind;

// ------------------------- labddu & wagddu ------------------------- //

wagdduadv=etaadv*wagddsadv+(1-etaadv)*waguadv;
wagddunam=etanam*wagddsnam+(1-etanam)*wagunam;
wagddussa=etassa*wagddsssa+(1-etassa)*wagussa;
wagddulac=etalac*wagddslac+(1-etalac)*wagulac;
wagddujap=etajap*wagddsjap+(1-etajap)*wagujap;
wagddueas=etaeas*wagddseas+(1-etaeas)*wagueas;
wagddumen=etamen*wagddsmen+(1-etamen)*wagumen;
wagddurus=etarus*wagddsrus+(1-etarus)*wagurus;
wagdduchi=etachi*wagddschi+(1-etachi)*waguchi;
wagdduind=etaind*wagddsind+(1-etaind)*waguind;


wagdduadv=(1-aaadv)*(1-alpha)*kadv^alpha*tfpadv^(1-alpha)*labdduadv^(sigma-1)
*(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^((1-alpha-sigma)/sigma);
wagddunam=(1-aanam)*(1-alpha)*knam^alpha*labddunam^(sigma-1)
*(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^((1-alpha-sigma)/sigma);
wagddussa=(1-aassa)*(1-alpha)*kssa^alpha*tfpssa^(1-alpha)*labddussa^(sigma-1)
*(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^((1-alpha-sigma)/sigma);
wagddulac=(1-aalac)*(1-alpha)*klac^alpha*tfplac^(1-alpha)*labddulac^(sigma-1)
*(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^((1-alpha-sigma)/sigma);
wagddujap=(1-aajap)*(1-alpha)*kjap^alpha*tfpjap^(1-alpha)*labddujap^(sigma-1)
*(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^((1-alpha-sigma)/sigma);
wagddueas=(1-aaeas)*(1-alpha)*keas^alpha*tfpeas^(1-alpha)*labddueas^(sigma-1)
*(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^((1-alpha-sigma)/sigma);
wagddumen=(1-aamen)*(1-alpha)*kmen^alpha*tfpmen^(1-alpha)*labddumen^(sigma-1)
*(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^((1-alpha-sigma)/sigma);
wagddurus=(1-aarus)*(1-alpha)*krus^alpha*tfprus^(1-alpha)*labddurus^(sigma-1)
*(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^((1-alpha-sigma)/sigma);
wagdduchi=(1-aachi)*(1-alpha)*kchi^alpha*tfpchi^(1-alpha)*labdduchi^(sigma-1)
*(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^((1-alpha-sigma)/sigma);
wagdduind=(1-aaind)*(1-alpha)*kind^alpha*tfpind^(1-alpha)*labdduind^(sigma-1)
*(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^((1-alpha-sigma)/sigma);

// --------------------------- unemeployement ----------------------- //

uradv=(labuadv-labdduadv)/labuadv;
urnam=(labunam-labddunam)/labunam;
urssa=(labussa-labddussa)/labussa;
urlac=(labulac-labddulac)/labulac;
urjap=(labujap-labddujap)/labujap;
ureas=(labueas-labddueas)/labueas;
urmen=(labumen-labddumen)/labumen;
urrus=(laburus-labddurus)/laburus;
urchi=(labuchi-labdduchi)/labuchi;
urind=(labuind-labdduind)/labuind;

// -------------------------------- ces ------------------------------- //

labadv=(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^(1/sigma);
labnam=(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^(1/sigma);
labssa=(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^(1/sigma);
lablac=(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^(1/sigma);
labjap=(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^(1/sigma);
labeas=(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^(1/sigma);
labmen=(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^(1/sigma);
labrus=(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^(1/sigma);
labchi=(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^(1/sigma);
labind=(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^(1/sigma);

// ----------------------------- wagdds ------------------------------- //

wagddsadv=aaadv*(1-alpha)*kadv^alpha*tfpadv^(1-alpha)*labsadv^(sigma-1)
*(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^((1-alpha-sigma)/sigma);
wagddsnam=aanam*(1-alpha)*knam^alpha*labsnam^(sigma-1)
*(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^((1-alpha-sigma)/sigma);
wagddsssa=aassa*(1-alpha)*kssa^alpha*tfpssa^(1-alpha)*labsssa^(sigma-1)
*(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^((1-alpha-sigma)/sigma);
wagddslac=aalac*(1-alpha)*klac^alpha*tfplac^(1-alpha)*labslac^(sigma-1)
*(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^((1-alpha-sigma)/sigma);
wagddsjap=aajap*(1-alpha)*kjap^alpha*tfpjap^(1-alpha)*labsjap^(sigma-1)
*(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^((1-alpha-sigma)/sigma);
wagddseas=aaeas*(1-alpha)*keas^alpha*tfpeas^(1-alpha)*labseas^(sigma-1)
*(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^((1-alpha-sigma)/sigma);
wagddsmen=aamen*(1-alpha)*kmen^alpha*tfpmen^(1-alpha)*labsmen^(sigma-1)
*(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^((1-alpha-sigma)/sigma);
wagddsrus=aarus*(1-alpha)*krus^alpha*tfprus^(1-alpha)*labsrus^(sigma-1)
*(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^((1-alpha-sigma)/sigma);
wagddschi=aachi*(1-alpha)*kchi^alpha*tfpchi^(1-alpha)*labschi^(sigma-1)
*(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^((1-alpha-sigma)/sigma);
wagddsind=aaind*(1-alpha)*kind^alpha*tfpind^(1-alpha)*labsind^(sigma-1)
*(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^((1-alpha-sigma)/sigma);

// ---------------------- wagddu and wagdds ------------------------- //

wagsadv=aaadv*(1-alpha)*kadv^alpha*tfpadv^(1-alpha)*labsadv^(sigma-1)
*(aaadv*labsadv^sigma+(1-aaadv)*labuadv^sigma)^((1-alpha-sigma)/sigma);
wagsnam=aanam*(1-alpha)*knam^alpha*labsnam^(sigma-1)
*(aanam*labsnam^sigma+(1-aanam)*labunam^sigma)^((1-alpha-sigma)/sigma);
wagsssa=aassa*(1-alpha)*kssa^alpha*tfpssa^(1-alpha)*labsssa^(sigma-1)
*(aassa*labsssa^sigma+(1-aassa)*labussa^sigma)^((1-alpha-sigma)/sigma);
wagslac=aalac*(1-alpha)*klac^alpha*tfplac^(1-alpha)*labslac^(sigma-1)
*(aalac*labslac^sigma+(1-aalac)*labulac^sigma)^((1-alpha-sigma)/sigma);
wagsjap=aajap*(1-alpha)*kjap^alpha*tfpjap^(1-alpha)*labsjap^(sigma-1)
*(aajap*labsjap^sigma+(1-aajap)*labujap^sigma)^((1-alpha-sigma)/sigma);
wagseas=aaeas*(1-alpha)*keas^alpha*tfpeas^(1-alpha)*labseas^(sigma-1)
*(aaeas*labseas^sigma+(1-aaeas)*labueas^sigma)^((1-alpha-sigma)/sigma);
wagsmen=aamen*(1-alpha)*kmen^alpha*tfpmen^(1-alpha)*labsmen^(sigma-1)
*(aamen*labsmen^sigma+(1-aamen)*labumen^sigma)^((1-alpha-sigma)/sigma);
wagsrus=aarus*(1-alpha)*krus^alpha*tfprus^(1-alpha)*labsrus^(sigma-1)
*(aarus*labsrus^sigma+(1-aarus)*laburus^sigma)^((1-alpha-sigma)/sigma);
wagschi=aachi*(1-alpha)*kchi^alpha*tfpchi^(1-alpha)*labschi^(sigma-1)
*(aachi*labschi^sigma+(1-aachi)*labuchi^sigma)^((1-alpha-sigma)/sigma);
wagsind=aaind*(1-alpha)*kind^alpha*tfpind^(1-alpha)*labsind^(sigma-1)
*(aaind*labsind^sigma+(1-aaind)*labuind^sigma)^((1-alpha-sigma)/sigma);

waguadv=(1-aaadv)*(1-alpha)*kadv^alpha*tfpadv^(1-alpha)*labuadv^(sigma-1)
*(aaadv*labsadv^sigma+(1-aaadv)*labuadv^sigma)^((1-alpha-sigma)/sigma);
wagunam=(1-aanam)*(1-alpha)*knam^alpha*labunam^(sigma-1)
*(aanam*labsnam^sigma+(1-aanam)*labunam^sigma)^((1-alpha-sigma)/sigma);
wagussa=(1-aassa)*(1-alpha)*kssa^alpha*tfpssa^(1-alpha)*labussa^(sigma-1)
*(aassa*labsssa^sigma+(1-aassa)*labussa^sigma)^((1-alpha-sigma)/sigma);
wagulac=(1-aalac)*(1-alpha)*klac^alpha*tfplac^(1-alpha)*labulac^(sigma-1)
*(aalac*labslac^sigma+(1-aalac)*labulac^sigma)^((1-alpha-sigma)/sigma);
wagujap=(1-aajap)*(1-alpha)*kjap^alpha*tfpjap^(1-alpha)*labujap^(sigma-1)
*(aajap*labsjap^sigma+(1-aajap)*labujap^sigma)^((1-alpha-sigma)/sigma);
wagueas=(1-aaeas)*(1-alpha)*keas^alpha*tfpeas^(1-alpha)*labueas^(sigma-1)
*(aaeas*labseas^sigma+(1-aaeas)*labueas^sigma)^((1-alpha-sigma)/sigma);
wagumen=(1-aamen)*(1-alpha)*kmen^alpha*tfpmen^(1-alpha)*labumen^(sigma-1)
*(aamen*labsmen^sigma+(1-aamen)*labumen^sigma)^((1-alpha-sigma)/sigma);
wagurus=(1-aarus)*(1-alpha)*krus^alpha*tfprus^(1-alpha)*laburus^(sigma-1)
*(aarus*labsrus^sigma+(1-aarus)*laburus^sigma)^((1-alpha-sigma)/sigma);
waguchi=(1-aachi)*(1-alpha)*kchi^alpha*tfpchi^(1-alpha)*labuchi^(sigma-1)
*(aachi*labschi^sigma+(1-aachi)*labuchi^sigma)^((1-alpha-sigma)/sigma);
waguind=(1-aaind)*(1-alpha)*kind^alpha*tfpind^(1-alpha)*labuind^(sigma-1)
*(aaind*labsind^sigma+(1-aaind)*labuind^sigma)^((1-alpha-sigma)/sigma);

//------------------------- ratiogdp -------------------------------- //

ratiogdpadv=gdpadv/gdpnam
/(1+p2adv/mmadv(-1)+p3adv/(mmadv(-1)*mmadv(-2))+p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))+p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdpssa=gdpssa/gdpnam
/(1+p2ssa/mmssa(-1)+p3ssa/(mmssa(-1)*mmssa(-2))+p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))+p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdplac=gdplac/gdpnam
/(1+p2lac/mmlac(-1)+p3lac/(mmlac(-1)*mmlac(-2))+p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))+p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdpjap=gdpjap/gdpnam
/(1+p2jap/mmjap(-1)+p3jap/(mmjap(-1)*mmjap(-2))+p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))+p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdprus=gdprus/gdpnam
/(1+p2rus/mmrus(-1)+p3rus/(mmrus(-1)*mmrus(-2))+p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))+p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdpmen=gdpmen/gdpnam
/(1+p2men/mmmen(-1)+p3men/(mmmen(-1)*mmmen(-2))+p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))+p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdpeas=gdpeas/gdpnam
/(1+p2eas/mmeas(-1)+p3eas/(mmeas(-1)*mmeas(-2))+p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))+p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdpchi=gdpchi/gdpnam
/(1+p2chi/mmchi(-1)+p3chi/(mmchi(-1)*mmchi(-2))+p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))+p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

ratiogdpind=gdpind/gdpnam
/(1+p2ind/mmind(-1)+p3ind/(mmind(-1)*mmind(-2))+p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))+p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
*(1+p2nam/mmnam(-1)+p3nam/(mmnam(-1)*mmnam(-2))+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))+p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));


// ----------------------- advanced countries ------------------------- //

xs8adv*(1+tcadv)=1/(p8adv*phiadv(-7))*intrate*zs8adv*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)+lams8adv*(1-tauadv)*wagddsadv+psiadv*tr8sadv*wagddsadv+(1-lams8adv)*bsadv;
xu8adv*(1+tcadv)=1/(p8adv*(1-phiadv(-7)))*intrate*zu8adv*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)+(1-uradv)*lamu8adv*(1-tauadv)*wagdduadv+psiadv*tr8uadv*wagdduadv+(1-lamu8adv)*buadv;

xs2adv(+1)*ggnam(+1)=beta*intrate(+1)*xs1adv*(1+tcadv)/(1+tcadv(+1));
xs3adv(+1)*ggnam(+1)=beta*intrate(+1)*xs2adv*(1+tcadv)/(1+tcadv(+1));
xs4adv(+1)*ggnam(+1)=beta*intrate(+1)*xs3adv*(1+tcadv)/(1+tcadv(+1));
xs5adv(+1)*ggnam(+1)=beta*intrate(+1)*xs4adv*(1+tcadv)/(1+tcadv(+1));
xs6adv(+1)*ggnam(+1)=beta*intrate(+1)*xs5adv*(1+tcadv)/(1+tcadv(+1));
xs7adv(+1)*ggnam(+1)=beta*intrate(+1)*xs6adv*(1+tcadv)/(1+tcadv(+1));
xs8adv(+1)*ggnam(+1)=beta*intrate(+1)*xs7adv*(1+tcadv)/(1+tcadv(+1));

xu2adv(+1)*ggnam(+1)=beta*intrate(+1)*xu1adv*(1+tcadv)/(1+tcadv(+1));
xu3adv(+1)*ggnam(+1)=beta*intrate(+1)*xu2adv*(1+tcadv)/(1+tcadv(+1));
xu4adv(+1)*ggnam(+1)=beta*intrate(+1)*xu3adv*(1+tcadv)/(1+tcadv(+1));
xu5adv(+1)*ggnam(+1)=beta*intrate(+1)*xu4adv*(1+tcadv)/(1+tcadv(+1));
xu6adv(+1)*ggnam(+1)=beta*intrate(+1)*xu5adv*(1+tcadv)/(1+tcadv(+1));
xu7adv(+1)*ggnam(+1)=beta*intrate(+1)*xu6adv*(1+tcadv)/(1+tcadv(+1));
xu8adv(+1)*ggnam(+1)=beta*intrate(+1)*xu7adv*(1+tcadv)/(1+tcadv(+1));

zs2adv=phiadv(-1)*((1-edusadv(-1))*(1-tauadv(-1))*lams1adv(-1)*wagddsadv(-1)+psiadv(-1)*tr1sadv(-1)*wagddsadv(-1)+(1-edusadv(-1))*(1-lams1adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs1adv(-1))/(ggnam*mmadv(-1));
zs3adv=intrate(-1)*zs2adv(-1)/(ggnam*mmadv(-1))+phiadv(-2)*p2adv(-1)*((1-tauadv(-1))*lams2adv(-1)*wagddsadv(-1)+psiadv(-1)*tr2sadv(-1)*wagddsadv(-1)+(1-lams2adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs2adv(-1))/(ggnam*mmadv(-1)*mmadv(-2));
zs4adv=intrate(-1)*zs3adv(-1)/(ggnam*mmadv(-1))+phiadv(-3)*p3adv(-1)*((1-tauadv(-1))*lams3adv(-1)*wagddsadv(-1)+psiadv(-1)*tr3sadv(-1)*wagddsadv(-1)+(1-lams3adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs3adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3));
zs5adv=intrate(-1)*zs4adv(-1)/(ggnam*mmadv(-1))+phiadv(-4)*p4adv(-1)*((1-tauadv(-1))*lams4adv(-1)*wagddsadv(-1)+psiadv(-1)*tr4sadv(-1)*wagddsadv(-1)+(1-lams4adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs4adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4));
zs6adv=intrate(-1)*zs5adv(-1)/(ggnam*mmadv(-1))+phiadv(-5)*p5adv(-1)*((1-tauadv(-1))*lams5adv(-1)*wagddsadv(-1)+psiadv(-1)*tr5sadv(-1)*wagddsadv(-1)+(1-lams5adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs5adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5));
zs7adv=intrate(-1)*zs6adv(-1)/(ggnam*mmadv(-1))+phiadv(-6)*p6adv(-1)*((1-tauadv(-1))*lams6adv(-1)*wagddsadv(-1)+psiadv(-1)*tr6sadv(-1)*wagddsadv(-1)+(1-lams6adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs6adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6));
zs8adv=intrate(-1)*zs7adv(-1)/(ggnam*mmadv(-1))+phiadv(-7)*p7adv(-1)*((1-tauadv(-1))*lams7adv(-1)*wagddsadv(-1)+psiadv(-1)*tr7sadv(-1)*wagddsadv(-1)+(1-lams7adv(-1))*bsadv(-1)-(1+tcadv(-1))*xs7adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7));

zu2adv=(1-phiadv(-1))*((1-eduuadv(-1))*(1-tauadv(-1))*lamu1adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr1uadv(-1)*wagdduadv(-1)+(1-eduuadv(-1))*(1-lamu1adv(-1))*buadv(-1)-(1+tcadv(-1))*xu1adv(-1))/(ggnam*mmadv(-1));
zu3adv=intrate(-1)*zu2adv(-1)/(ggnam*mmadv(-1))+(1-phiadv(-2))*p2adv(-1)*((1-tauadv(-1))*lamu2adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr2uadv(-1)*wagdduadv(-1)+(1-lamu2adv(-1))*buadv(-1)-(1+tcadv(-1))*xu2adv(-1))/(ggnam*mmadv(-1)*mmadv(-2));
zu4adv=intrate(-1)*zu3adv(-1)/(ggnam*mmadv(-1))+(1-phiadv(-3))*p3adv(-1)*((1-tauadv(-1))*lamu3adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr3uadv(-1)*wagdduadv(-1)+(1-lamu3adv(-1))*buadv(-1)-(1+tcadv(-1))*xu3adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3));
zu5adv=intrate(-1)*zu4adv(-1)/(ggnam*mmadv(-1))+(1-phiadv(-4))*p4adv(-1)*((1-tauadv(-1))*lamu4adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr4uadv(-1)*wagdduadv(-1)+(1-lamu4adv(-1))*buadv(-1)-(1+tcadv(-1))*xu4adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4));
zu6adv=intrate(-1)*zu5adv(-1)/(ggnam*mmadv(-1))+(1-phiadv(-5))*p5adv(-1)*((1-tauadv(-1))*lamu5adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr5uadv(-1)*wagdduadv(-1)+(1-lamu5adv(-1))*buadv(-1)-(1+tcadv(-1))*xu5adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5));
zu7adv=intrate(-1)*zu6adv(-1)/(ggnam*mmadv(-1))+(1-phiadv(-6))*p6adv(-1)*((1-tauadv(-1))*lamu6adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr6uadv(-1)*wagdduadv(-1)+(1-lamu6adv(-1))*buadv(-1)-(1+tcadv(-1))*xu6adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6));
zu8adv=intrate(-1)*zu7adv(-1)/(ggnam*mmadv(-1))+(1-phiadv(-7))*p7adv(-1)*((1-tauadv(-1))*lamu7adv(-1)*(1-uradv)*wagdduadv(-1)+psiadv(-1)*tr7uadv(-1)*wagdduadv(-1)+(1-lamu7adv(-1))*buadv(-1)-(1+tcadv(-1))*xu7adv(-1))/(ggnam*mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7));

weaadv=zs2adv+zs3adv+zs4adv+zs5adv+zs6adv+zs7adv+zs8adv
+zu2adv+zu3adv+zu4adv+zu5adv+zu6adv+zu7adv+zu8adv;

gdpadv=kadv^alpha*(tfpadv*labadv)^(1-alpha);

govadv=tauadv*(wagddsadv*labsadv+wagdduadv*labdduadv)
+tcadv*(phiadv*xs1adv
+xs2adv*phiadv(-1)*p2adv/mmadv(-1)
+xs3adv*phiadv(-2)*p3adv/(mmadv(-1)*mmadv(-2))
+xs4adv*phiadv(-3)*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+xs5adv*phiadv(-4)*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+xs6adv*phiadv(-5)*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+xs7adv*phiadv(-6)*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+xs8adv*phiadv(-7)*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
+tcadv*((1-phiadv)*xu1adv
+xu2adv*(1-phiadv(-1))*p2adv/mmadv(-1)
+xu3adv*(1-phiadv(-2))*p3adv/(mmadv(-1)*mmadv(-2))
+xu4adv*(1-phiadv(-3))*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+xu5adv*(1-phiadv(-4))*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+xu6adv*(1-phiadv(-5))*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+xu7adv*(1-phiadv(-6))*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+xu8adv*(1-phiadv(-7))*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
-bsadv*((1-lams1adv)*phiadv*(1-edusadv)
+(1-lams2adv)*phiadv(-1)*p2adv/mmadv(-1)
+(1-lams3adv)*phiadv(-2)*p3adv/(mmadv(-1)*mmadv(-2))
+(1-lams4adv)*phiadv(-3)*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+(1-lams5adv)*phiadv(-4)*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+(1-lams6adv)*phiadv(-5)*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+(1-lams7adv)*phiadv(-6)*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+(1-lams8adv)*phiadv(-7)*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
-buadv*((1-lamu1adv)*(1-phiadv)*(1-eduuadv)
+(1-lamu2adv)*(1-phiadv(-1))*p2adv/mmadv(-1)
+(1-lamu3adv)*(1-phiadv(-2))*p3adv/(mmadv(-1)*mmadv(-2))
+(1-lamu4adv)*(1-phiadv(-3))*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+(1-lamu5adv)*(1-phiadv(-4))*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+(1-lamu6adv)*(1-phiadv(-5))*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+(1-lamu7adv)*(1-phiadv(-6))*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+(1-lamu8adv)*(1-phiadv(-7))*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
-psiadv*wagddsadv*(tr1sadv*phiadv
+tr2sadv*phiadv(-1)*p2adv/mmadv(-1)
+tr3sadv*phiadv(-2)*p3adv/(mmadv(-1)*mmadv(-2))
+tr4sadv*phiadv(-3)*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+tr5sadv*phiadv(-4)*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+tr6sadv*phiadv(-5)*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+tr7sadv*phiadv(-6)*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+tr8sadv*phiadv(-7)*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
-psiadv*wagdduadv*(tr1uadv*(1-phiadv)
+tr2uadv*(1-phiadv(-1))*p2adv/mmadv(-1)
+tr3uadv*(1-phiadv(-2))*p3adv/(mmadv(-1)*mmadv(-2))
+tr4uadv*(1-phiadv(-3))*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+tr5uadv*(1-phiadv(-4))*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+tr6uadv*(1-phiadv(-5))*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+tr7uadv*(1-phiadv(-6))*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+tr8uadv*(1-phiadv(-7))*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
-conspubgdpadv;

// -------------------------------- nam ------------------------------- //

xs8nam*(1+tcnam)=1/(p8nam*phinam(-7))*intrate*zs8nam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)+lams8nam*(1-taunam)*wagddsnam+psinam*tr8snam*wagddsnam+(1-lams8nam)*bsnam;
xu8nam*(1+tcnam)=1/(p8nam*(1-phinam(-7)))*intrate*zu8nam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)+(1-urnam)*lamu8nam*(1-taunam)*wagddunam+psinam*tr8unam*wagddunam+(1-lamu8nam)*bunam;

xs2nam(+1)*ggnam(+1)=beta*intrate(+1)*xs1nam*(1+tcnam)/(1+tcnam(+1));
xs3nam(+1)*ggnam(+1)=beta*intrate(+1)*xs2nam*(1+tcnam)/(1+tcnam(+1));
xs4nam(+1)*ggnam(+1)=beta*intrate(+1)*xs3nam*(1+tcnam)/(1+tcnam(+1));
xs5nam(+1)*ggnam(+1)=beta*intrate(+1)*xs4nam*(1+tcnam)/(1+tcnam(+1));
xs6nam(+1)*ggnam(+1)=beta*intrate(+1)*xs5nam*(1+tcnam)/(1+tcnam(+1));
xs7nam(+1)*ggnam(+1)=beta*intrate(+1)*xs6nam*(1+tcnam)/(1+tcnam(+1));
xs8nam(+1)*ggnam(+1)=beta*intrate(+1)*xs7nam*(1+tcnam)/(1+tcnam(+1));

xu2nam(+1)*ggnam(+1)=beta*intrate(+1)*xu1nam*(1+tcnam)/(1+tcnam(+1));
xu3nam(+1)*ggnam(+1)=beta*intrate(+1)*xu2nam*(1+tcnam)/(1+tcnam(+1));
xu4nam(+1)*ggnam(+1)=beta*intrate(+1)*xu3nam*(1+tcnam)/(1+tcnam(+1));
xu5nam(+1)*ggnam(+1)=beta*intrate(+1)*xu4nam*(1+tcnam)/(1+tcnam(+1));
xu6nam(+1)*ggnam(+1)=beta*intrate(+1)*xu5nam*(1+tcnam)/(1+tcnam(+1));
xu7nam(+1)*ggnam(+1)=beta*intrate(+1)*xu6nam*(1+tcnam)/(1+tcnam(+1));
xu8nam(+1)*ggnam(+1)=beta*intrate(+1)*xu7nam*(1+tcnam)/(1+tcnam(+1));

zs2nam=phinam(-1)*((1-edusnam(-1))*(1-taunam(-1))*lams1nam(-1)*wagddsnam(-1)+psinam(-1)*tr1snam(-1)*wagddsnam(-1)+(1-edusnam(-1))*(1-lams1nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs1nam(-1))/(ggnam*mmnam(-1));
zs3nam=intrate(-1)*zs2nam(-1)/(ggnam*mmnam(-1))+phinam(-2)*p2nam(-1)*((1-taunam(-1))*lams2nam(-1)*wagddsnam(-1)+psinam(-1)*tr2snam(-1)*wagddsnam(-1)+(1-lams2nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs2nam(-1))/(ggnam*mmnam(-1)*mmnam(-2));
zs4nam=intrate(-1)*zs3nam(-1)/(ggnam*mmnam(-1))+phinam(-3)*p3nam(-1)*((1-taunam(-1))*lams3nam(-1)*wagddsnam(-1)+psinam(-1)*tr3snam(-1)*wagddsnam(-1)+(1-lams3nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs3nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3));
zs5nam=intrate(-1)*zs4nam(-1)/(ggnam*mmnam(-1))+phinam(-4)*p4nam(-1)*((1-taunam(-1))*lams4nam(-1)*wagddsnam(-1)+psinam(-1)*tr4snam(-1)*wagddsnam(-1)+(1-lams4nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs4nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4));
zs6nam=intrate(-1)*zs5nam(-1)/(ggnam*mmnam(-1))+phinam(-5)*p5nam(-1)*((1-taunam(-1))*lams5nam(-1)*wagddsnam(-1)+psinam(-1)*tr5snam(-1)*wagddsnam(-1)+(1-lams5nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs5nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5));
zs7nam=intrate(-1)*zs6nam(-1)/(ggnam*mmnam(-1))+phinam(-6)*p6nam(-1)*((1-taunam(-1))*lams6nam(-1)*wagddsnam(-1)+psinam(-1)*tr6snam(-1)*wagddsnam(-1)+(1-lams6nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs6nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6));
zs8nam=intrate(-1)*zs7nam(-1)/(ggnam*mmnam(-1))+phinam(-7)*p7nam(-1)*((1-taunam(-1))*lams7nam(-1)*wagddsnam(-1)+psinam(-1)*tr7snam(-1)*wagddsnam(-1)+(1-lams7nam(-1))*bsnam(-1)-(1+tcnam(-1))*xs7nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7));

zu2nam=(1-phinam(-1))*((1-eduunam(-1))*(1-taunam(-1))*lamu1nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr1unam(-1)*wagddunam(-1)+(1-eduunam(-1))*(1-lamu1nam(-1))*bunam(-1)-(1+tcnam(-1))*xu1nam(-1))/(ggnam*mmnam(-1));
zu3nam=intrate(-1)*zu2nam(-1)/(ggnam*mmnam(-1))+(1-phinam(-2))*p2nam(-1)*((1-taunam(-1))*lamu2nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr2unam(-1)*wagddunam(-1)+(1-lamu2nam(-1))*bunam(-1)-(1+tcnam(-1))*xu2nam(-1))/(ggnam*mmnam(-1)*mmnam(-2));
zu4nam=intrate(-1)*zu3nam(-1)/(ggnam*mmnam(-1))+(1-phinam(-3))*p3nam(-1)*((1-taunam(-1))*lamu3nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr3unam(-1)*wagddunam(-1)+(1-lamu3nam(-1))*bunam(-1)-(1+tcnam(-1))*xu3nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3));
zu5nam=intrate(-1)*zu4nam(-1)/(ggnam*mmnam(-1))+(1-phinam(-4))*p4nam(-1)*((1-taunam(-1))*lamu4nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr4unam(-1)*wagddunam(-1)+(1-lamu4nam(-1))*bunam(-1)-(1+tcnam(-1))*xu4nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4));
zu6nam=intrate(-1)*zu5nam(-1)/(ggnam*mmnam(-1))+(1-phinam(-5))*p5nam(-1)*((1-taunam(-1))*lamu5nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr5unam(-1)*wagddunam(-1)+(1-lamu5nam(-1))*bunam(-1)-(1+tcnam(-1))*xu5nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5));
zu7nam=intrate(-1)*zu6nam(-1)/(ggnam*mmnam(-1))+(1-phinam(-6))*p6nam(-1)*((1-taunam(-1))*lamu6nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr6unam(-1)*wagddunam(-1)+(1-lamu6nam(-1))*bunam(-1)-(1+tcnam(-1))*xu6nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6));
zu8nam=intrate(-1)*zu7nam(-1)/(ggnam*mmnam(-1))+(1-phinam(-7))*p7nam(-1)*((1-taunam(-1))*lamu7nam(-1)*(1-urnam)*wagddunam(-1)+psinam(-1)*tr7unam(-1)*wagddunam(-1)+(1-lamu7nam(-1))*bunam(-1)-(1+tcnam(-1))*xu7nam(-1))/(ggnam*mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7));

weanam=zs2nam+zs3nam+zs4nam+zs5nam+zs6nam+zs7nam+zs8nam
+zu2nam+zu3nam+zu4nam+zu5nam+zu6nam+zu7nam+zu8nam;

gdpnam=knam^alpha*labnam^(1-alpha);

govnam=taunam*(wagddsnam*labsnam+wagddunam*labddunam)
+tcnam*(phinam*xs1nam
+xs2nam*phinam(-1)*p2nam/mmnam(-1)
+xs3nam*phinam(-2)*p3nam/(mmnam(-1)*mmnam(-2))
+xs4nam*phinam(-3)*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+xs5nam*phinam(-4)*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+xs6nam*phinam(-5)*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+xs7nam*phinam(-6)*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+xs8nam*phinam(-7)*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
+tcnam*((1-phinam)*xu1nam
+xu2nam*(1-phinam(-1))*p2nam/mmnam(-1)
+xu3nam*(1-phinam(-2))*p3nam/(mmnam(-1)*mmnam(-2))
+xu4nam*(1-phinam(-3))*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+xu5nam*(1-phinam(-4))*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+xu6nam*(1-phinam(-5))*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+xu7nam*(1-phinam(-6))*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+xu8nam*(1-phinam(-7))*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
-bsnam*((1-lams1nam)*phinam*(1-edusnam)
+(1-lams2nam)*phinam(-1)*p2nam/mmnam(-1)
+(1-lams3nam)*phinam(-2)*p3nam/(mmnam(-1)*mmnam(-2))
+(1-lams4nam)*phinam(-3)*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+(1-lams5nam)*phinam(-4)*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+(1-lams6nam)*phinam(-5)*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+(1-lams7nam)*phinam(-6)*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+(1-lams8nam)*phinam(-7)*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
-bunam*((1-lamu1nam)*(1-phinam)*(1-eduunam)
+(1-lamu2nam)*(1-phinam(-1))*p2nam/mmnam(-1)
+(1-lamu3nam)*(1-phinam(-2))*p3nam/(mmnam(-1)*mmnam(-2))
+(1-lamu4nam)*(1-phinam(-3))*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+(1-lamu5nam)*(1-phinam(-4))*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+(1-lamu6nam)*(1-phinam(-5))*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+(1-lamu7nam)*(1-phinam(-6))*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+(1-lamu8nam)*(1-phinam(-7))*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
-psinam*wagddsnam*(tr1snam*phinam
+tr2snam*phinam(-1)*p2nam/mmnam(-1)
+tr3snam*phinam(-2)*p3nam/(mmnam(-1)*mmnam(-2))
+tr4snam*phinam(-3)*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+tr5snam*phinam(-4)*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+tr6snam*phinam(-5)*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+tr7snam*phinam(-6)*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+tr8snam*phinam(-7)*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
-psinam*wagddunam*(tr1unam*(1-phinam)
+tr2unam*(1-phinam(-1))*p2nam/mmnam(-1)
+tr3unam*(1-phinam(-2))*p3nam/(mmnam(-1)*mmnam(-2))
+tr4unam*(1-phinam(-3))*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+tr5unam*(1-phinam(-4))*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+tr6unam*(1-phinam(-5))*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+tr7unam*(1-phinam(-6))*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+tr8unam*(1-phinam(-7))*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
-conspubgdpnam;

// ------------------------------ japan ------------------------------- //

xs8jap*(1+tcjap)=1/(p8jap*phijap(-7))*intrate*zs8jap*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)+lams8jap*(1-taujap)*wagddsjap+psijap*tr8sjap*wagddsjap+(1-lams8jap)*bsjap;
xu8jap*(1+tcjap)=1/(p8jap*(1-phijap(-7)))*intrate*zu8jap*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)+(1-urjap)*lamu8jap*(1-taujap)*wagddujap+psijap*tr8ujap*wagddujap+(1-lamu8jap)*bujap;

xs2jap(+1)*ggnam(+1)=beta*intrate(+1)*xs1jap*(1+tcjap)/(1+tcjap(+1));
xs3jap(+1)*ggnam(+1)=beta*intrate(+1)*xs2jap*(1+tcjap)/(1+tcjap(+1));
xs4jap(+1)*ggnam(+1)=beta*intrate(+1)*xs3jap*(1+tcjap)/(1+tcjap(+1));
xs5jap(+1)*ggnam(+1)=beta*intrate(+1)*xs4jap*(1+tcjap)/(1+tcjap(+1));
xs6jap(+1)*ggnam(+1)=beta*intrate(+1)*xs5jap*(1+tcjap)/(1+tcjap(+1));
xs7jap(+1)*ggnam(+1)=beta*intrate(+1)*xs6jap*(1+tcjap)/(1+tcjap(+1));
xs8jap(+1)*ggnam(+1)=beta*intrate(+1)*xs7jap*(1+tcjap)/(1+tcjap(+1));

xu2jap(+1)*ggnam(+1)=beta*intrate(+1)*xu1jap*(1+tcjap)/(1+tcjap(+1));
xu3jap(+1)*ggnam(+1)=beta*intrate(+1)*xu2jap*(1+tcjap)/(1+tcjap(+1));
xu4jap(+1)*ggnam(+1)=beta*intrate(+1)*xu3jap*(1+tcjap)/(1+tcjap(+1));
xu5jap(+1)*ggnam(+1)=beta*intrate(+1)*xu4jap*(1+tcjap)/(1+tcjap(+1));
xu6jap(+1)*ggnam(+1)=beta*intrate(+1)*xu5jap*(1+tcjap)/(1+tcjap(+1));
xu7jap(+1)*ggnam(+1)=beta*intrate(+1)*xu6jap*(1+tcjap)/(1+tcjap(+1));
xu8jap(+1)*ggnam(+1)=beta*intrate(+1)*xu7jap*(1+tcjap)/(1+tcjap(+1));

zs2jap=phijap(-1)*((1-edusjap(-1))*(1-taujap(-1))*lams1jap(-1)*wagddsjap(-1)+psijap(-1)*tr1sjap(-1)*wagddsjap(-1)+(1-edusjap(-1))*(1-lams1jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs1jap(-1))/(ggnam*mmjap(-1));
zs3jap=intrate(-1)*zs2jap(-1)/(ggnam*mmjap(-1))+phijap(-2)*p2jap(-1)*((1-taujap(-1))*lams2jap(-1)*wagddsjap(-1)+psijap(-1)*tr2sjap(-1)*wagddsjap(-1)+(1-lams2jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs2jap(-1))/(ggnam*mmjap(-1)*mmjap(-2));
zs4jap=intrate(-1)*zs3jap(-1)/(ggnam*mmjap(-1))+phijap(-3)*p3jap(-1)*((1-taujap(-1))*lams3jap(-1)*wagddsjap(-1)+psijap(-1)*tr3sjap(-1)*wagddsjap(-1)+(1-lams3jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs3jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3));
zs5jap=intrate(-1)*zs4jap(-1)/(ggnam*mmjap(-1))+phijap(-4)*p4jap(-1)*((1-taujap(-1))*lams4jap(-1)*wagddsjap(-1)+psijap(-1)*tr4sjap(-1)*wagddsjap(-1)+(1-lams4jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs4jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4));
zs6jap=intrate(-1)*zs5jap(-1)/(ggnam*mmjap(-1))+phijap(-5)*p5jap(-1)*((1-taujap(-1))*lams5jap(-1)*wagddsjap(-1)+psijap(-1)*tr5sjap(-1)*wagddsjap(-1)+(1-lams5jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs5jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5));
zs7jap=intrate(-1)*zs6jap(-1)/(ggnam*mmjap(-1))+phijap(-6)*p6jap(-1)*((1-taujap(-1))*lams6jap(-1)*wagddsjap(-1)+psijap(-1)*tr6sjap(-1)*wagddsjap(-1)+(1-lams6jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs6jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6));
zs8jap=intrate(-1)*zs7jap(-1)/(ggnam*mmjap(-1))+phijap(-7)*p7jap(-1)*((1-taujap(-1))*lams7jap(-1)*wagddsjap(-1)+psijap(-1)*tr7sjap(-1)*wagddsjap(-1)+(1-lams7jap(-1))*bsjap(-1)-(1+tcjap(-1))*xs7jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7));

zu2jap=(1-phijap(-1))*((1-eduujap(-1))*(1-taujap(-1))*lamu1jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr1ujap(-1)*wagddujap(-1)+(1-eduujap(-1))*(1-lamu1jap(-1))*bujap(-1)-(1+tcjap(-1))*xu1jap(-1))/(ggnam*mmjap(-1));
zu3jap=intrate(-1)*zu2jap(-1)/(ggnam*mmjap(-1))+(1-phijap(-2))*p2jap(-1)*((1-taujap(-1))*lamu2jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr2ujap(-1)*wagddujap(-1)+(1-lamu2jap(-1))*bujap(-1)-(1+tcjap(-1))*xu2jap(-1))/(ggnam*mmjap(-1)*mmjap(-2));
zu4jap=intrate(-1)*zu3jap(-1)/(ggnam*mmjap(-1))+(1-phijap(-3))*p3jap(-1)*((1-taujap(-1))*lamu3jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr3ujap(-1)*wagddujap(-1)+(1-lamu3jap(-1))*bujap(-1)-(1+tcjap(-1))*xu3jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3));
zu5jap=intrate(-1)*zu4jap(-1)/(ggnam*mmjap(-1))+(1-phijap(-4))*p4jap(-1)*((1-taujap(-1))*lamu4jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr4ujap(-1)*wagddujap(-1)+(1-lamu4jap(-1))*bujap(-1)-(1+tcjap(-1))*xu4jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4));
zu6jap=intrate(-1)*zu5jap(-1)/(ggnam*mmjap(-1))+(1-phijap(-5))*p5jap(-1)*((1-taujap(-1))*lamu5jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr5ujap(-1)*wagddujap(-1)+(1-lamu5jap(-1))*bujap(-1)-(1+tcjap(-1))*xu5jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5));
zu7jap=intrate(-1)*zu6jap(-1)/(ggnam*mmjap(-1))+(1-phijap(-6))*p6jap(-1)*((1-taujap(-1))*lamu6jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr6ujap(-1)*wagddujap(-1)+(1-lamu6jap(-1))*bujap(-1)-(1+tcjap(-1))*xu6jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6));
zu8jap=intrate(-1)*zu7jap(-1)/(ggnam*mmjap(-1))+(1-phijap(-7))*p7jap(-1)*((1-taujap(-1))*lamu7jap(-1)*(1-urjap(-1))*wagddujap(-1)+psijap(-1)*tr7ujap(-1)*wagddujap(-1)+(1-lamu7jap(-1))*bujap(-1)-(1+tcjap(-1))*xu7jap(-1))/(ggnam*mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7));

weajap=zs2jap+zs3jap+zs4jap+zs5jap+zs6jap+zs7jap+zs8jap
+zu2jap+zu3jap+zu4jap+zu5jap+zu6jap+zu7jap+zu8jap;

gdpjap=kjap^alpha*(tfpjap*labjap)^(1-alpha);

govjap=taujap*(wagddsjap*labsjap+wagddujap*labddujap)
+tcjap*(phijap*xs1jap
+xs2jap*phijap(-1)*p2jap/mmjap(-1)
+xs3jap*phijap(-2)*p3jap/(mmjap(-1)*mmjap(-2))
+xs4jap*phijap(-3)*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+xs5jap*phijap(-4)*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+xs6jap*phijap(-5)*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+xs7jap*phijap(-6)*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+xs8jap*phijap(-7)*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
+tcjap*((1-phijap)*xu1jap
+xu2jap*(1-phijap(-1))*p2jap/mmjap(-1)
+xu3jap*(1-phijap(-2))*p3jap/(mmjap(-1)*mmjap(-2))
+xu4jap*(1-phijap(-3))*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+xu5jap*(1-phijap(-4))*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+xu6jap*(1-phijap(-5))*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+xu7jap*(1-phijap(-6))*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+xu8jap*(1-phijap(-7))*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
-bsjap*((1-lams1jap)*phijap*(1-edusjap)
+(1-lams2jap)*phijap(-1)*p2jap/mmjap(-1)
+(1-lams3jap)*phijap(-2)*p3jap/(mmjap(-1)*mmjap(-2))
+(1-lams4jap)*phijap(-3)*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+(1-lams5jap)*phijap(-4)*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+(1-lams6jap)*phijap(-5)*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+(1-lams7jap)*phijap(-6)*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+(1-lams8jap)*phijap(-7)*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
-bujap*((1-lamu1jap)*(1-phijap)*(1-eduujap)
+(1-lamu2jap)*(1-phijap(-1))*p2jap/mmjap(-1)
+(1-lamu3jap)*(1-phijap(-2))*p3jap/(mmjap(-1)*mmjap(-2))
+(1-lamu4jap)*(1-phijap(-3))*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+(1-lamu5jap)*(1-phijap(-4))*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+(1-lamu6jap)*(1-phijap(-5))*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+(1-lamu7jap)*(1-phijap(-6))*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+(1-lamu8jap)*(1-phijap(-7))*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
-psijap*wagddsjap*(tr1sjap*phijap
+tr2sjap*phijap(-1)*p2jap/mmjap(-1)
+tr3sjap*phijap(-2)*p3jap/(mmjap(-1)*mmjap(-2))
+tr4sjap*phijap(-3)*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+tr5sjap*phijap(-4)*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+tr6sjap*phijap(-5)*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+tr7sjap*phijap(-6)*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+tr8sjap*phijap(-7)*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
-psijap*wagddujap*(tr1ujap*(1-phijap)
+tr2ujap*(1-phijap(-1))*p2jap/mmjap(-1)
+tr3ujap*(1-phijap(-2))*p3jap/(mmjap(-1)*mmjap(-2))
+tr4ujap*(1-phijap(-3))*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+tr5ujap*(1-phijap(-4))*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+tr6ujap*(1-phijap(-5))*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+tr7ujap*(1-phijap(-6))*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+tr8ujap*(1-phijap(-7))*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
-conspubgdpjap;

// -------------------------------- ssa ------------------------------- //

xs8ssa*(1+tcssa)=1/(p8ssa*phissa(-7))*intrate*zs8ssa*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)+lams8ssa*(1-taussa)*wagddsssa+psissa*tr8sssa*wagddsssa+(1-lams8ssa)*bsssa;
xu8ssa*(1+tcssa)=1/(p8ssa*(1-phissa(-7)))*intrate*zu8ssa*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)+(1-urssa)*lamu8ssa*(1-taussa)*wagddussa+psissa*tr8ussa*wagddussa+(1-lamu8ssa)*bussa;

xs2ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs1ssa*(1+tcssa)/(1+tcssa(+1));
xs3ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs2ssa*(1+tcssa)/(1+tcssa(+1));
xs4ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs3ssa*(1+tcssa)/(1+tcssa(+1));
xs5ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs4ssa*(1+tcssa)/(1+tcssa(+1));
xs6ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs5ssa*(1+tcssa)/(1+tcssa(+1));
xs7ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs6ssa*(1+tcssa)/(1+tcssa(+1));
xs8ssa(+1)*ggnam(+1)=beta*intrate(+1)*xs7ssa*(1+tcssa)/(1+tcssa(+1));

xu2ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu1ssa*(1+tcssa)/(1+tcssa(+1));
xu3ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu2ssa*(1+tcssa)/(1+tcssa(+1));
xu4ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu3ssa*(1+tcssa)/(1+tcssa(+1));
xu5ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu4ssa*(1+tcssa)/(1+tcssa(+1));
xu6ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu5ssa*(1+tcssa)/(1+tcssa(+1));
xu7ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu6ssa*(1+tcssa)/(1+tcssa(+1));
xu8ssa(+1)*ggnam(+1)=beta*intrate(+1)*xu7ssa*(1+tcssa)/(1+tcssa(+1));

zs2ssa=phissa(-1)*((1-edusssa(-1))*(1-taussa(-1))*lams1ssa(-1)*wagddsssa(-1)+psissa(-1)*tr1sssa(-1)*wagddsssa(-1)+(1-edusssa(-1))*(1-lams1ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs1ssa(-1))/(ggnam*mmssa(-1));
zs3ssa=intrate(-1)*zs2ssa(-1)/(ggnam*mmssa(-1))+phissa(-2)*p2ssa(-1)*((1-taussa(-1))*lams2ssa(-1)*wagddsssa(-1)+psissa(-1)*tr2sssa(-1)*wagddsssa(-1)+(1-lams2ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs2ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2));
zs4ssa=intrate(-1)*zs3ssa(-1)/(ggnam*mmssa(-1))+phissa(-3)*p3ssa(-1)*((1-taussa(-1))*lams3ssa(-1)*wagddsssa(-1)+psissa(-1)*tr3sssa(-1)*wagddsssa(-1)+(1-lams3ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs3ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3));
zs5ssa=intrate(-1)*zs4ssa(-1)/(ggnam*mmssa(-1))+phissa(-4)*p4ssa(-1)*((1-taussa(-1))*lams4ssa(-1)*wagddsssa(-1)+psissa(-1)*tr4sssa(-1)*wagddsssa(-1)+(1-lams4ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs4ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4));
zs6ssa=intrate(-1)*zs5ssa(-1)/(ggnam*mmssa(-1))+phissa(-5)*p5ssa(-1)*((1-taussa(-1))*lams5ssa(-1)*wagddsssa(-1)+psissa(-1)*tr5sssa(-1)*wagddsssa(-1)+(1-lams5ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs5ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5));
zs7ssa=intrate(-1)*zs6ssa(-1)/(ggnam*mmssa(-1))+phissa(-6)*p6ssa(-1)*((1-taussa(-1))*lams6ssa(-1)*wagddsssa(-1)+psissa(-1)*tr6sssa(-1)*wagddsssa(-1)+(1-lams6ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs6ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6));
zs8ssa=intrate(-1)*zs7ssa(-1)/(ggnam*mmssa(-1))+phissa(-7)*p7ssa(-1)*((1-taussa(-1))*lams7ssa(-1)*wagddsssa(-1)+psissa(-1)*tr7sssa(-1)*wagddsssa(-1)+(1-lams7ssa(-1))*bsssa(-1)-(1+tcssa(-1))*xs7ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7));

zu2ssa=(1-phissa(-1))*((1-eduussa(-1))*(1-taussa(-1))*lamu1ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr1ussa(-1)*wagddussa(-1)+(1-eduussa(-1))*(1-lamu1ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu1ssa(-1))/(ggnam*mmssa(-1));
zu3ssa=intrate(-1)*zu2ssa(-1)/(ggnam*mmssa(-1))+(1-phissa(-2))*p2ssa(-1)*((1-taussa(-1))*lamu2ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr2ussa(-1)*wagddussa(-1)+(1-lamu2ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu2ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2));
zu4ssa=intrate(-1)*zu3ssa(-1)/(ggnam*mmssa(-1))+(1-phissa(-3))*p3ssa(-1)*((1-taussa(-1))*lamu3ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr3ussa(-1)*wagddussa(-1)+(1-lamu3ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu3ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3));
zu5ssa=intrate(-1)*zu4ssa(-1)/(ggnam*mmssa(-1))+(1-phissa(-4))*p4ssa(-1)*((1-taussa(-1))*lamu4ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr4ussa(-1)*wagddussa(-1)+(1-lamu4ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu4ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4));
zu6ssa=intrate(-1)*zu5ssa(-1)/(ggnam*mmssa(-1))+(1-phissa(-5))*p5ssa(-1)*((1-taussa(-1))*lamu5ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr5ussa(-1)*wagddussa(-1)+(1-lamu5ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu5ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5));
zu7ssa=intrate(-1)*zu6ssa(-1)/(ggnam*mmssa(-1))+(1-phissa(-6))*p6ssa(-1)*((1-taussa(-1))*lamu6ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr6ussa(-1)*wagddussa(-1)+(1-lamu6ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu6ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6));
zu8ssa=intrate(-1)*zu7ssa(-1)/(ggnam*mmssa(-1))+(1-phissa(-7))*p7ssa(-1)*((1-taussa(-1))*lamu7ssa(-1)*(1-urssa(-1))*wagddussa(-1)+psissa(-1)*tr7ussa(-1)*wagddussa(-1)+(1-lamu7ssa(-1))*bussa(-1)-(1+tcssa(-1))*xu7ssa(-1))/(ggnam*mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7));

weassa=zs2ssa+zs3ssa+zs4ssa+zs5ssa+zs6ssa+zs7ssa+zs8ssa
+zu2ssa+zu3ssa+zu4ssa+zu5ssa+zu6ssa+zu7ssa+zu8ssa;

gdpssa=kssa^alpha*(tfpssa*labssa)^(1-alpha);

govssa=taussa*(wagddsssa*labsssa+wagddussa*labddussa)
+pissa*intrate*kssa
+tcssa*(phissa*xs1ssa
+xs2ssa*phissa(-1)*p2ssa/mmssa(-1)
+xs3ssa*phissa(-2)*p3ssa/(mmssa(-1)*mmssa(-2))
+xs4ssa*phissa(-3)*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+xs5ssa*phissa(-4)*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+xs6ssa*phissa(-5)*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+xs7ssa*phissa(-6)*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+xs8ssa*phissa(-7)*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
+tcssa*((1-phissa)*xu1ssa
+xu2ssa*(1-phissa(-1))*p2ssa/mmssa(-1)
+xu3ssa*(1-phissa(-2))*p3ssa/(mmssa(-1)*mmssa(-2))
+xu4ssa*(1-phissa(-3))*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+xu5ssa*(1-phissa(-4))*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+xu6ssa*(1-phissa(-5))*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+xu7ssa*(1-phissa(-6))*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+xu8ssa*(1-phissa(-7))*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
-bsssa*((1-lams1ssa)*phissa*(1-edusssa)
+(1-lams2ssa)*phissa(-1)*p2ssa/mmssa(-1)
+(1-lams3ssa)*phissa(-2)*p3ssa/(mmssa(-1)*mmssa(-2))
+(1-lams4ssa)*phissa(-3)*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+(1-lams5ssa)*phissa(-4)*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+(1-lams6ssa)*phissa(-5)*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+(1-lams7ssa)*phissa(-6)*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+(1-lams8ssa)*phissa(-7)*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
-bussa*((1-lamu1ssa)*(1-phissa)*(1-eduussa)
+(1-lamu2ssa)*(1-phissa(-1))*p2ssa/mmssa(-1)
+(1-lamu3ssa)*(1-phissa(-2))*p3ssa/(mmssa(-1)*mmssa(-2))
+(1-lamu4ssa)*(1-phissa(-3))*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+(1-lamu5ssa)*(1-phissa(-4))*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+(1-lamu6ssa)*(1-phissa(-5))*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+(1-lamu7ssa)*(1-phissa(-6))*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+(1-lamu8ssa)*(1-phissa(-7))*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
-psissa*wagddsssa*(tr1sssa*phissa
+tr2sssa*phissa(-1)*p2ssa/mmssa(-1)
+tr3sssa*phissa(-2)*p3ssa/(mmssa(-1)*mmssa(-2))
+tr4sssa*phissa(-3)*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+tr5sssa*phissa(-4)*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+tr6sssa*phissa(-5)*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+tr7sssa*phissa(-6)*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+tr8sssa*phissa(-7)*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
-psissa*wagddussa*(tr1ussa*(1-phissa)
+tr2ussa*(1-phissa(-1))*p2ssa/mmssa(-1)
+tr3ussa*(1-phissa(-2))*p3ssa/(mmssa(-1)*mmssa(-2))
+tr4ussa*(1-phissa(-3))*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+tr5ussa*(1-phissa(-4))*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+tr6ussa*(1-phissa(-5))*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+tr7ussa*(1-phissa(-6))*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+tr8ussa*(1-phissa(-7))*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
-conspubgdpssa
-pissa*intrate*kssa;

// ------------------------------- lac ------------------------------- //

xs8lac*(1+tclac)=1/(p8lac*philac(-7))*intrate*zs8lac*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)+lams8lac*(1-taulac)*wagddslac+psilac*tr8slac*wagddslac+(1-lams8lac)*bslac;
xu8lac*(1+tclac)=1/(p8lac*(1-philac(-7)))*intrate*zu8lac*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)+(1-urlac)*lamu8lac*(1-taulac)*wagddulac+psilac*tr8ulac*wagddulac+(1-lamu8lac)*bulac;

xs2lac(+1)*ggnam(+1)=beta*intrate(+1)*xs1lac*(1+tclac)/(1+tclac(+1));
xs3lac(+1)*ggnam(+1)=beta*intrate(+1)*xs2lac*(1+tclac)/(1+tclac(+1));
xs4lac(+1)*ggnam(+1)=beta*intrate(+1)*xs3lac*(1+tclac)/(1+tclac(+1));
xs5lac(+1)*ggnam(+1)=beta*intrate(+1)*xs4lac*(1+tclac)/(1+tclac(+1));
xs6lac(+1)*ggnam(+1)=beta*intrate(+1)*xs5lac*(1+tclac)/(1+tclac(+1));
xs7lac(+1)*ggnam(+1)=beta*intrate(+1)*xs6lac*(1+tclac)/(1+tclac(+1));
xs8lac(+1)*ggnam(+1)=beta*intrate(+1)*xs7lac*(1+tclac)/(1+tclac(+1));

xu2lac(+1)*ggnam(+1)=beta*intrate(+1)*xu1lac*(1+tclac)/(1+tclac(+1));
xu3lac(+1)*ggnam(+1)=beta*intrate(+1)*xu2lac*(1+tclac)/(1+tclac(+1));
xu4lac(+1)*ggnam(+1)=beta*intrate(+1)*xu3lac*(1+tclac)/(1+tclac(+1));
xu5lac(+1)*ggnam(+1)=beta*intrate(+1)*xu4lac*(1+tclac)/(1+tclac(+1));
xu6lac(+1)*ggnam(+1)=beta*intrate(+1)*xu5lac*(1+tclac)/(1+tclac(+1));
xu7lac(+1)*ggnam(+1)=beta*intrate(+1)*xu6lac*(1+tclac)/(1+tclac(+1));
xu8lac(+1)*ggnam(+1)=beta*intrate(+1)*xu7lac*(1+tclac)/(1+tclac(+1));

zs2lac=philac(-1)*((1-eduslac(-1))*(1-taulac(-1))*lams1lac(-1)*wagddslac(-1)+psilac(-1)*tr1slac(-1)*wagddslac(-1)+(1-eduslac(-1))*(1-lams1lac(-1))*bslac(-1)-(1+tclac(-1))*xs1lac(-1))/(ggnam*mmlac(-1));
zs3lac=intrate(-1)*zs2lac(-1)/(ggnam*mmlac(-1))+philac(-2)*p2lac(-1)*((1-taulac(-1))*lams2lac(-1)*wagddslac(-1)+psilac(-1)*tr2slac(-1)*wagddslac(-1)+(1-lams2lac(-1))*bslac(-1)-(1+tclac(-1))*xs2lac(-1))/(ggnam*mmlac(-1)*mmlac(-2));
zs4lac=intrate(-1)*zs3lac(-1)/(ggnam*mmlac(-1))+philac(-3)*p3lac(-1)*((1-taulac(-1))*lams3lac(-1)*wagddslac(-1)+psilac(-1)*tr3slac(-1)*wagddslac(-1)+(1-lams3lac(-1))*bslac(-1)-(1+tclac(-1))*xs3lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3));
zs5lac=intrate(-1)*zs4lac(-1)/(ggnam*mmlac(-1))+philac(-4)*p4lac(-1)*((1-taulac(-1))*lams4lac(-1)*wagddslac(-1)+psilac(-1)*tr4slac(-1)*wagddslac(-1)+(1-lams4lac(-1))*bslac(-1)-(1+tclac(-1))*xs4lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4));
zs6lac=intrate(-1)*zs5lac(-1)/(ggnam*mmlac(-1))+philac(-5)*p5lac(-1)*((1-taulac(-1))*lams5lac(-1)*wagddslac(-1)+psilac(-1)*tr5slac(-1)*wagddslac(-1)+(1-lams5lac(-1))*bslac(-1)-(1+tclac(-1))*xs5lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5));
zs7lac=intrate(-1)*zs6lac(-1)/(ggnam*mmlac(-1))+philac(-6)*p6lac(-1)*((1-taulac(-1))*lams6lac(-1)*wagddslac(-1)+psilac(-1)*tr6slac(-1)*wagddslac(-1)+(1-lams6lac(-1))*bslac(-1)-(1+tclac(-1))*xs6lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6));
zs8lac=intrate(-1)*zs7lac(-1)/(ggnam*mmlac(-1))+philac(-7)*p7lac(-1)*((1-taulac(-1))*lams7lac(-1)*wagddslac(-1)+psilac(-1)*tr7slac(-1)*wagddslac(-1)+(1-lams7lac(-1))*bslac(-1)-(1+tclac(-1))*xs7lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7));

zu2lac=(1-philac(-1))*((1-eduulac(-1))*(1-taulac(-1))*lamu1lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr1ulac(-1)*wagddulac(-1)+(1-eduulac(-1))*(1-lamu1lac(-1))*bulac(-1)-(1+tclac(-1))*xu1lac(-1))/(ggnam*mmlac(-1));
zu3lac=intrate(-1)*zu2lac(-1)/(ggnam*mmlac(-1))+(1-philac(-2))*p2lac(-1)*((1-taulac(-1))*lamu2lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr2ulac(-1)*wagddulac(-1)+(1-lamu2lac(-1))*bulac(-1)-(1+tclac(-1))*xu2lac(-1))/(ggnam*mmlac(-1)*mmlac(-2));
zu4lac=intrate(-1)*zu3lac(-1)/(ggnam*mmlac(-1))+(1-philac(-3))*p3lac(-1)*((1-taulac(-1))*lamu3lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr3ulac(-1)*wagddulac(-1)+(1-lamu3lac(-1))*bulac(-1)-(1+tclac(-1))*xu3lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3));
zu5lac=intrate(-1)*zu4lac(-1)/(ggnam*mmlac(-1))+(1-philac(-4))*p4lac(-1)*((1-taulac(-1))*lamu4lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr4ulac(-1)*wagddulac(-1)+(1-lamu4lac(-1))*bulac(-1)-(1+tclac(-1))*xu4lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4));
zu6lac=intrate(-1)*zu5lac(-1)/(ggnam*mmlac(-1))+(1-philac(-5))*p5lac(-1)*((1-taulac(-1))*lamu5lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr5ulac(-1)*wagddulac(-1)+(1-lamu5lac(-1))*bulac(-1)-(1+tclac(-1))*xu5lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5));
zu7lac=intrate(-1)*zu6lac(-1)/(ggnam*mmlac(-1))+(1-philac(-6))*p6lac(-1)*((1-taulac(-1))*lamu6lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr6ulac(-1)*wagddulac(-1)+(1-lamu6lac(-1))*bulac(-1)-(1+tclac(-1))*xu6lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6));
zu8lac=intrate(-1)*zu7lac(-1)/(ggnam*mmlac(-1))+(1-philac(-7))*p7lac(-1)*((1-taulac(-1))*lamu7lac(-1)*(1-urlac(-1))*wagddulac(-1)+psilac(-1)*tr7ulac(-1)*wagddulac(-1)+(1-lamu7lac(-1))*bulac(-1)-(1+tclac(-1))*xu7lac(-1))/(ggnam*mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7));

wealac=zs2lac+zs3lac+zs4lac+zs5lac+zs6lac+zs7lac+zs8lac
+zu2lac+zu3lac+zu4lac+zu5lac+zu6lac+zu7lac+zu8lac;

gdplac=klac^alpha*(tfplac*lablac)^(1-alpha);

govlac=taulac*(wagddslac*labslac+wagddulac*labddulac)
+pilac*intrate*klac
+tclac*(philac*xs1lac
+xs2lac*philac(-1)*p2lac/mmlac(-1)
+xs3lac*philac(-2)*p3lac/(mmlac(-1)*mmlac(-2))
+xs4lac*philac(-3)*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+xs5lac*philac(-4)*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+xs6lac*philac(-5)*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+xs7lac*philac(-6)*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+xs8lac*philac(-7)*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
+tclac*((1-philac)*xu1lac
+xu2lac*(1-philac(-1))*p2lac/mmlac(-1)
+xu3lac*(1-philac(-2))*p3lac/(mmlac(-1)*mmlac(-2))
+xu4lac*(1-philac(-3))*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+xu5lac*(1-philac(-4))*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+xu6lac*(1-philac(-5))*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+xu7lac*(1-philac(-6))*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+xu8lac*(1-philac(-7))*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
-bslac*((1-lams1lac)*philac*(1-eduslac)
+(1-lams2lac)*philac(-1)*p2lac/mmlac(-1)
+(1-lams3lac)*philac(-2)*p3lac/(mmlac(-1)*mmlac(-2))
+(1-lams4lac)*philac(-3)*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+(1-lams5lac)*philac(-4)*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+(1-lams6lac)*philac(-5)*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+(1-lams7lac)*philac(-6)*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+(1-lams8lac)*philac(-7)*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
-bulac*((1-lamu1lac)*(1-philac)*(1-eduulac)
+(1-lamu2lac)*(1-philac(-1))*p2lac/mmlac(-1)
+(1-lamu3lac)*(1-philac(-2))*p3lac/(mmlac(-1)*mmlac(-2))
+(1-lamu4lac)*(1-philac(-3))*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+(1-lamu5lac)*(1-philac(-4))*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+(1-lamu6lac)*(1-philac(-5))*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+(1-lamu7lac)*(1-philac(-6))*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+(1-lamu8lac)*(1-philac(-7))*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
-psilac*wagddslac*(tr1slac*philac
+tr2slac*philac(-1)*p2lac/mmlac(-1)
+tr3slac*philac(-2)*p3lac/(mmlac(-1)*mmlac(-2))
+tr4slac*philac(-3)*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+tr5slac*philac(-4)*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+tr6slac*philac(-5)*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+tr7slac*philac(-6)*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+tr8slac*philac(-7)*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
-psilac*wagddulac*(tr1ulac*(1-philac)
+tr2ulac*(1-philac(-1))*p2lac/mmlac(-1)
+tr3ulac*(1-philac(-2))*p3lac/(mmlac(-1)*mmlac(-2))
+tr4ulac*(1-philac(-3))*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+tr5ulac*(1-philac(-4))*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+tr6ulac*(1-philac(-5))*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+tr7ulac*(1-philac(-6))*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+tr8ulac*(1-philac(-7))*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
-conspubgdplac
-pilac*intrate*klac;

// ---------------------------- rusland ------------------------------- //

xs8rus*(1+tcrus)=1/(p8rus*phirus(-7))*intrate*zs8rus*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)+lams8rus*(1-taurus)*wagddsrus+psirus*tr8srus*wagddsrus+(1-lams8rus)*bsrus;
xu8rus*(1+tcrus)=1/(p8rus*(1-phirus(-7)))*intrate*zu8rus*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)+(1-urrus)*lamu8rus*(1-taurus)*wagddurus+psirus*tr8urus*wagddurus+(1-lamu8rus)*burus;

xs2rus(+1)*ggnam(+1)=beta*intrate(+1)*xs1rus*(1+tcrus)/(1+tcrus(+1));
xs3rus(+1)*ggnam(+1)=beta*intrate(+1)*xs2rus*(1+tcrus)/(1+tcrus(+1));
xs4rus(+1)*ggnam(+1)=beta*intrate(+1)*xs3rus*(1+tcrus)/(1+tcrus(+1));
xs5rus(+1)*ggnam(+1)=beta*intrate(+1)*xs4rus*(1+tcrus)/(1+tcrus(+1));
xs6rus(+1)*ggnam(+1)=beta*intrate(+1)*xs5rus*(1+tcrus)/(1+tcrus(+1));
xs7rus(+1)*ggnam(+1)=beta*intrate(+1)*xs6rus*(1+tcrus)/(1+tcrus(+1));
xs8rus(+1)*ggnam(+1)=beta*intrate(+1)*xs7rus*(1+tcrus)/(1+tcrus(+1));

xu2rus(+1)*ggnam(+1)=beta*intrate(+1)*xu1rus*(1+tcrus)/(1+tcrus(+1));
xu3rus(+1)*ggnam(+1)=beta*intrate(+1)*xu2rus*(1+tcrus)/(1+tcrus(+1));
xu4rus(+1)*ggnam(+1)=beta*intrate(+1)*xu3rus*(1+tcrus)/(1+tcrus(+1));
xu5rus(+1)*ggnam(+1)=beta*intrate(+1)*xu4rus*(1+tcrus)/(1+tcrus(+1));
xu6rus(+1)*ggnam(+1)=beta*intrate(+1)*xu5rus*(1+tcrus)/(1+tcrus(+1));
xu7rus(+1)*ggnam(+1)=beta*intrate(+1)*xu6rus*(1+tcrus)/(1+tcrus(+1));
xu8rus(+1)*ggnam(+1)=beta*intrate(+1)*xu7rus*(1+tcrus)/(1+tcrus(+1));

zs2rus=phirus(-1)*((1-edusrus(-1))*(1-taurus(-1))*lams1rus(-1)*wagddsrus(-1)+psirus(-1)*tr1srus(-1)*wagddsrus(-1)+(1-edusrus(-1))*(1-lams1rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs1rus(-1))/(ggnam*mmrus(-1));
zs3rus=intrate(-1)*zs2rus(-1)/(ggnam*mmrus(-1))+phirus(-2)*p2rus(-1)*((1-taurus(-1))*lams2rus(-1)*wagddsrus(-1)+psirus(-1)*tr2srus(-1)*wagddsrus(-1)+(1-lams2rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs2rus(-1))/(ggnam*mmrus(-1)*mmrus(-2));
zs4rus=intrate(-1)*zs3rus(-1)/(ggnam*mmrus(-1))+phirus(-3)*p3rus(-1)*((1-taurus(-1))*lams3rus(-1)*wagddsrus(-1)+psirus(-1)*tr3srus(-1)*wagddsrus(-1)+(1-lams3rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs3rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3));
zs5rus=intrate(-1)*zs4rus(-1)/(ggnam*mmrus(-1))+phirus(-4)*p4rus(-1)*((1-taurus(-1))*lams4rus(-1)*wagddsrus(-1)+psirus(-1)*tr4srus(-1)*wagddsrus(-1)+(1-lams4rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs4rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4));
zs6rus=intrate(-1)*zs5rus(-1)/(ggnam*mmrus(-1))+phirus(-5)*p5rus(-1)*((1-taurus(-1))*lams5rus(-1)*wagddsrus(-1)+psirus(-1)*tr5srus(-1)*wagddsrus(-1)+(1-lams5rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs5rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5));
zs7rus=intrate(-1)*zs6rus(-1)/(ggnam*mmrus(-1))+phirus(-6)*p6rus(-1)*((1-taurus(-1))*lams6rus(-1)*wagddsrus(-1)+psirus(-1)*tr6srus(-1)*wagddsrus(-1)+(1-lams6rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs6rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6));
zs8rus=intrate(-1)*zs7rus(-1)/(ggnam*mmrus(-1))+phirus(-7)*p7rus(-1)*((1-taurus(-1))*lams7rus(-1)*wagddsrus(-1)+psirus(-1)*tr7srus(-1)*wagddsrus(-1)+(1-lams7rus(-1))*bsrus(-1)-(1+tcrus(-1))*xs7rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7));

zu2rus=(1-phirus(-1))*((1-eduurus(-1))*(1-taurus(-1))*lamu1rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr1urus(-1)*wagddurus(-1)+(1-eduurus(-1))*(1-lamu1rus(-1))*burus(-1)-(1+tcrus(-1))*xu1rus(-1))/(ggnam*mmrus(-1));
zu3rus=intrate(-1)*zu2rus(-1)/(ggnam*mmrus(-1))+(1-phirus(-2))*p2rus(-1)*((1-taurus(-1))*lamu2rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr2urus(-1)*wagddurus(-1)+(1-lamu2rus(-1))*burus(-1)-(1+tcrus(-1))*xu2rus(-1))/(ggnam*mmrus(-1)*mmrus(-2));
zu4rus=intrate(-1)*zu3rus(-1)/(ggnam*mmrus(-1))+(1-phirus(-3))*p3rus(-1)*((1-taurus(-1))*lamu3rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr3urus(-1)*wagddurus(-1)+(1-lamu3rus(-1))*burus(-1)-(1+tcrus(-1))*xu3rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3));
zu5rus=intrate(-1)*zu4rus(-1)/(ggnam*mmrus(-1))+(1-phirus(-4))*p4rus(-1)*((1-taurus(-1))*lamu4rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr4urus(-1)*wagddurus(-1)+(1-lamu4rus(-1))*burus(-1)-(1+tcrus(-1))*xu4rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4));
zu6rus=intrate(-1)*zu5rus(-1)/(ggnam*mmrus(-1))+(1-phirus(-5))*p5rus(-1)*((1-taurus(-1))*lamu5rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr5urus(-1)*wagddurus(-1)+(1-lamu5rus(-1))*burus(-1)-(1+tcrus(-1))*xu5rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5));
zu7rus=intrate(-1)*zu6rus(-1)/(ggnam*mmrus(-1))+(1-phirus(-6))*p6rus(-1)*((1-taurus(-1))*lamu6rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr6urus(-1)*wagddurus(-1)+(1-lamu6rus(-1))*burus(-1)-(1+tcrus(-1))*xu6rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6));
zu8rus=intrate(-1)*zu7rus(-1)/(ggnam*mmrus(-1))+(1-phirus(-7))*p7rus(-1)*((1-taurus(-1))*lamu7rus(-1)*(1-urrus)*wagddurus(-1)+psirus(-1)*tr7urus(-1)*wagddurus(-1)+(1-lamu7rus(-1))*burus(-1)-(1+tcrus(-1))*xu7rus(-1))/(ggnam*mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7));

wearus=zs2rus+zs3rus+zs4rus+zs5rus+zs6rus+zs7rus+zs8rus
+zu2rus+zu3rus+zu4rus+zu5rus+zu6rus+zu7rus+zu8rus;

gdprus=krus^alpha*(tfprus*labrus)^(1-alpha);

govrus=taurus*(wagddsrus*labsrus+wagddurus*labddurus)
+pirus*intrate*krus
+tcrus*(phirus*xs1rus
+xs2rus*phirus(-1)*p2rus/mmrus(-1)
+xs3rus*phirus(-2)*p3rus/(mmrus(-1)*mmrus(-2))
+xs4rus*phirus(-3)*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+xs5rus*phirus(-4)*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+xs6rus*phirus(-5)*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+xs7rus*phirus(-6)*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+xs8rus*phirus(-7)*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
+tcrus*((1-phirus)*xu1rus
+xu2rus*(1-phirus(-1))*p2rus/mmrus(-1)
+xu3rus*(1-phirus(-2))*p3rus/(mmrus(-1)*mmrus(-2))
+xu4rus*(1-phirus(-3))*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+xu5rus*(1-phirus(-4))*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+xu6rus*(1-phirus(-5))*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+xu7rus*(1-phirus(-6))*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+xu8rus*(1-phirus(-7))*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
-bsrus*((1-lams1rus)*phirus*(1-edusrus)
+(1-lams2rus)*phirus(-1)*p2rus/mmrus(-1)
+(1-lams3rus)*phirus(-2)*p3rus/(mmrus(-1)*mmrus(-2))
+(1-lams4rus)*phirus(-3)*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+(1-lams5rus)*phirus(-4)*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+(1-lams6rus)*phirus(-5)*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+(1-lams7rus)*phirus(-6)*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+(1-lams8rus)*phirus(-7)*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
-burus*((1-lamu1rus)*(1-phirus)*(1-eduurus)
+(1-lamu2rus)*(1-phirus(-1))*p2rus/mmrus(-1)
+(1-lamu3rus)*(1-phirus(-2))*p3rus/(mmrus(-1)*mmrus(-2))
+(1-lamu4rus)*(1-phirus(-3))*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+(1-lamu5rus)*(1-phirus(-4))*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+(1-lamu6rus)*(1-phirus(-5))*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+(1-lamu7rus)*(1-phirus(-6))*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+(1-lamu8rus)*(1-phirus(-7))*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
-psirus*wagddsrus*(tr1srus*phirus
+tr2srus*phirus(-1)*p2rus/mmrus(-1)
+tr3srus*phirus(-2)*p3rus/(mmrus(-1)*mmrus(-2))
+tr4srus*phirus(-3)*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+tr5srus*phirus(-4)*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+tr6srus*phirus(-5)*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+tr7srus*phirus(-6)*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+tr8srus*phirus(-7)*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
-psirus*wagddurus*(tr1urus*(1-phirus)
+tr2urus*(1-phirus(-1))*p2rus/mmrus(-1)
+tr3urus*(1-phirus(-2))*p3rus/(mmrus(-1)*mmrus(-2))
+tr4urus*(1-phirus(-3))*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+tr5urus*(1-phirus(-4))*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+tr6urus*(1-phirus(-5))*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+tr7urus*(1-phirus(-6))*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+tr8urus*(1-phirus(-7))*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
-conspubgdprus
-pirus*intrate*krus;

// ------------------------------ mena -------------------------------- //

xs8men*(1+tcmen)=1/(p8men*phimen(-7))*intrate*zs8men*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)+lams8men*(1-taumen)*wagddsmen+psimen*tr8smen*wagddsmen+(1-lams8men)*bsmen;
xu8men*(1+tcmen)=1/(p8men*(1-phimen(-7)))*intrate*zu8men*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)+(1-urmen)*lamu8men*(1-taumen)*wagddumen+psimen*tr8umen*wagddumen+(1-lamu8men)*bumen;

xs2men(+1)*ggnam(+1)=beta*intrate(+1)*xs1men*(1+tcmen)/(1+tcmen(+1));
xs3men(+1)*ggnam(+1)=beta*intrate(+1)*xs2men*(1+tcmen)/(1+tcmen(+1));
xs4men(+1)*ggnam(+1)=beta*intrate(+1)*xs3men*(1+tcmen)/(1+tcmen(+1));
xs5men(+1)*ggnam(+1)=beta*intrate(+1)*xs4men*(1+tcmen)/(1+tcmen(+1));
xs6men(+1)*ggnam(+1)=beta*intrate(+1)*xs5men*(1+tcmen)/(1+tcmen(+1));
xs7men(+1)*ggnam(+1)=beta*intrate(+1)*xs6men*(1+tcmen)/(1+tcmen(+1));
xs8men(+1)*ggnam(+1)=beta*intrate(+1)*xs7men*(1+tcmen)/(1+tcmen(+1));

xu2men(+1)*ggnam(+1)=beta*intrate(+1)*xu1men*(1+tcmen)/(1+tcmen(+1));
xu3men(+1)*ggnam(+1)=beta*intrate(+1)*xu2men*(1+tcmen)/(1+tcmen(+1));
xu4men(+1)*ggnam(+1)=beta*intrate(+1)*xu3men*(1+tcmen)/(1+tcmen(+1));
xu5men(+1)*ggnam(+1)=beta*intrate(+1)*xu4men*(1+tcmen)/(1+tcmen(+1));
xu6men(+1)*ggnam(+1)=beta*intrate(+1)*xu5men*(1+tcmen)/(1+tcmen(+1));
xu7men(+1)*ggnam(+1)=beta*intrate(+1)*xu6men*(1+tcmen)/(1+tcmen(+1));
xu8men(+1)*ggnam(+1)=beta*intrate(+1)*xu7men*(1+tcmen)/(1+tcmen(+1));

zs2men=phimen(-1)*((1-edusmen(-1))*(1-taumen(-1))*lams1men(-1)*wagddsmen(-1)+psimen(-1)*tr1smen(-1)*wagddsmen(-1)+(1-edusmen(-1))*(1-lams1men(-1))*bsmen(-1)-(1+tcmen(-1))*xs1men(-1))/(ggnam*mmmen(-1));
zs3men=intrate(-1)*zs2men(-1)/(ggnam*mmmen(-1))+phimen(-2)*p2men(-1)*((1-taumen(-1))*lams2men(-1)*wagddsmen(-1)+psimen(-1)*tr2smen(-1)*wagddsmen(-1)+(1-lams2men(-1))*bsmen(-1)-(1+tcmen(-1))*xs2men(-1))/(ggnam*mmmen(-1)*mmmen(-2));
zs4men=intrate(-1)*zs3men(-1)/(ggnam*mmmen(-1))+phimen(-3)*p3men(-1)*((1-taumen(-1))*lams3men(-1)*wagddsmen(-1)+psimen(-1)*tr3smen(-1)*wagddsmen(-1)+(1-lams3men(-1))*bsmen(-1)-(1+tcmen(-1))*xs3men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3));
zs5men=intrate(-1)*zs4men(-1)/(ggnam*mmmen(-1))+phimen(-4)*p4men(-1)*((1-taumen(-1))*lams4men(-1)*wagddsmen(-1)+psimen(-1)*tr4smen(-1)*wagddsmen(-1)+(1-lams4men(-1))*bsmen(-1)-(1+tcmen(-1))*xs4men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4));
zs6men=intrate(-1)*zs5men(-1)/(ggnam*mmmen(-1))+phimen(-5)*p5men(-1)*((1-taumen(-1))*lams5men(-1)*wagddsmen(-1)+psimen(-1)*tr5smen(-1)*wagddsmen(-1)+(1-lams5men(-1))*bsmen(-1)-(1+tcmen(-1))*xs5men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5));
zs7men=intrate(-1)*zs6men(-1)/(ggnam*mmmen(-1))+phimen(-6)*p6men(-1)*((1-taumen(-1))*lams6men(-1)*wagddsmen(-1)+psimen(-1)*tr6smen(-1)*wagddsmen(-1)+(1-lams6men(-1))*bsmen(-1)-(1+tcmen(-1))*xs6men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6));
zs8men=intrate(-1)*zs7men(-1)/(ggnam*mmmen(-1))+phimen(-7)*p7men(-1)*((1-taumen(-1))*lams7men(-1)*wagddsmen(-1)+psimen(-1)*tr7smen(-1)*wagddsmen(-1)+(1-lams7men(-1))*bsmen(-1)-(1+tcmen(-1))*xs7men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7));

zu2men=(1-phimen(-1))*((1-eduumen(-1))*(1-taumen(-1))*lamu1men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr1umen(-1)*wagddumen(-1)+(1-eduumen(-1))*(1-lamu1men(-1))*bumen(-1)-(1+tcmen(-1))*xu1men(-1))/(ggnam*mmmen(-1));
zu3men=intrate(-1)*zu2men(-1)/(ggnam*mmmen(-1))+(1-phimen(-2))*p2men(-1)*((1-taumen(-1))*lamu2men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr2umen(-1)*wagddumen(-1)+(1-lamu2men(-1))*bumen(-1)-(1+tcmen(-1))*xu2men(-1))/(ggnam*mmmen(-1)*mmmen(-2));
zu4men=intrate(-1)*zu3men(-1)/(ggnam*mmmen(-1))+(1-phimen(-3))*p3men(-1)*((1-taumen(-1))*lamu3men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr3umen(-1)*wagddumen(-1)+(1-lamu3men(-1))*bumen(-1)-(1+tcmen(-1))*xu3men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3));
zu5men=intrate(-1)*zu4men(-1)/(ggnam*mmmen(-1))+(1-phimen(-4))*p4men(-1)*((1-taumen(-1))*lamu4men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr4umen(-1)*wagddumen(-1)+(1-lamu4men(-1))*bumen(-1)-(1+tcmen(-1))*xu4men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4));
zu6men=intrate(-1)*zu5men(-1)/(ggnam*mmmen(-1))+(1-phimen(-5))*p5men(-1)*((1-taumen(-1))*lamu5men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr5umen(-1)*wagddumen(-1)+(1-lamu5men(-1))*bumen(-1)-(1+tcmen(-1))*xu5men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5));
zu7men=intrate(-1)*zu6men(-1)/(ggnam*mmmen(-1))+(1-phimen(-6))*p6men(-1)*((1-taumen(-1))*lamu6men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr6umen(-1)*wagddumen(-1)+(1-lamu6men(-1))*bumen(-1)-(1+tcmen(-1))*xu6men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6));
zu8men=intrate(-1)*zu7men(-1)/(ggnam*mmmen(-1))+(1-phimen(-7))*p7men(-1)*((1-taumen(-1))*lamu7men(-1)*(1-urmen(-1))*wagddumen(-1)+psimen(-1)*tr7umen(-1)*wagddumen(-1)+(1-lamu7men(-1))*bumen(-1)-(1+tcmen(-1))*xu7men(-1))/(ggnam*mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7));

weamen=zs2men+zs3men+zs4men+zs5men+zs6men+zs7men+zs8men
+zu2men+zu3men+zu4men+zu5men+zu6men+zu7men+zu8men;

gdpmen=kmen^alpha*(tfpmen*labmen)^(1-alpha);

govmen=taumen*(wagddsmen*labsmen+wagddumen*labddumen)
+pimen*intrate*kmen
+tcmen*(phimen*xs1men
+xs2men*phimen(-1)*p2men/mmmen(-1)
+xs3men*phimen(-2)*p3men/(mmmen(-1)*mmmen(-2))
+xs4men*phimen(-3)*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+xs5men*phimen(-4)*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+xs6men*phimen(-5)*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+xs7men*phimen(-6)*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+xs8men*phimen(-7)*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
+tcmen*((1-phimen)*xu1men
+xu2men*(1-phimen(-1))*p2men/mmmen(-1)
+xu3men*(1-phimen(-2))*p3men/(mmmen(-1)*mmmen(-2))
+xu4men*(1-phimen(-3))*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+xu5men*(1-phimen(-4))*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+xu6men*(1-phimen(-5))*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+xu7men*(1-phimen(-6))*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+xu8men*(1-phimen(-7))*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
-bsmen*((1-lams1men)*phimen*(1-edusmen)
+(1-lams2men)*phimen(-1)*p2men/mmmen(-1)
+(1-lams3men)*phimen(-2)*p3men/(mmmen(-1)*mmmen(-2))
+(1-lams4men)*phimen(-3)*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+(1-lams5men)*phimen(-4)*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+(1-lams6men)*phimen(-5)*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+(1-lams7men)*phimen(-6)*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+(1-lams8men)*phimen(-7)*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
-bumen*((1-lamu1men)*(1-phimen)*(1-eduumen)
+(1-lamu2men)*(1-phimen(-1))*p2men/mmmen(-1)
+(1-lamu3men)*(1-phimen(-2))*p3men/(mmmen(-1)*mmmen(-2))
+(1-lamu4men)*(1-phimen(-3))*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+(1-lamu5men)*(1-phimen(-4))*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+(1-lamu6men)*(1-phimen(-5))*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+(1-lamu7men)*(1-phimen(-6))*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+(1-lamu8men)*(1-phimen(-7))*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
-psimen*wagddsmen*(tr1smen*phimen
+tr2smen*phimen(-1)*p2men/mmmen(-1)
+tr3smen*phimen(-2)*p3men/(mmmen(-1)*mmmen(-2))
+tr4smen*phimen(-3)*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+tr5smen*phimen(-4)*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+tr6smen*phimen(-5)*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+tr7smen*phimen(-6)*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+tr8smen*phimen(-7)*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
-psimen*wagddumen*(tr1umen*(1-phimen)
+tr2umen*(1-phimen(-1))*p2men/mmmen(-1)
+tr3umen*(1-phimen(-2))*p3men/(mmmen(-1)*mmmen(-2))
+tr4umen*(1-phimen(-3))*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+tr5umen*(1-phimen(-4))*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+tr6umen*(1-phimen(-5))*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+tr7umen*(1-phimen(-6))*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+tr8umen*(1-phimen(-7))*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
-conspubgdpmen
-pimen*intrate*kmen;

// ------------------------ eastern europe ---------------------------- //

xs8eas*(1+tceas)=1/(p8eas*phieas(-7))*intrate*zs8eas*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)+lams8eas*(1-taueas)*wagddseas+psieas*tr8seas*wagddseas+(1-lams8eas)*bseas;
xu8eas*(1+tceas)=1/(p8eas*(1-phieas(-7)))*intrate*zu8eas*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)+(1-ureas)*lamu8eas*(1-taueas)*wagddueas+psieas*tr8ueas*wagddueas+(1-lamu8eas)*bueas;

xs2eas(+1)*ggnam(+1)=beta*intrate(+1)*xs1eas*(1+tceas)/(1+tceas(+1));
xs3eas(+1)*ggnam(+1)=beta*intrate(+1)*xs2eas*(1+tceas)/(1+tceas(+1));
xs4eas(+1)*ggnam(+1)=beta*intrate(+1)*xs3eas*(1+tceas)/(1+tceas(+1));
xs5eas(+1)*ggnam(+1)=beta*intrate(+1)*xs4eas*(1+tceas)/(1+tceas(+1));
xs6eas(+1)*ggnam(+1)=beta*intrate(+1)*xs5eas*(1+tceas)/(1+tceas(+1));
xs7eas(+1)*ggnam(+1)=beta*intrate(+1)*xs6eas*(1+tceas)/(1+tceas(+1));
xs8eas(+1)*ggnam(+1)=beta*intrate(+1)*xs7eas*(1+tceas)/(1+tceas(+1));

xu2eas(+1)*ggnam(+1)=beta*intrate(+1)*xu1eas*(1+tceas)/(1+tceas(+1));
xu3eas(+1)*ggnam(+1)=beta*intrate(+1)*xu2eas*(1+tceas)/(1+tceas(+1));
xu4eas(+1)*ggnam(+1)=beta*intrate(+1)*xu3eas*(1+tceas)/(1+tceas(+1));
xu5eas(+1)*ggnam(+1)=beta*intrate(+1)*xu4eas*(1+tceas)/(1+tceas(+1));
xu6eas(+1)*ggnam(+1)=beta*intrate(+1)*xu5eas*(1+tceas)/(1+tceas(+1));
xu7eas(+1)*ggnam(+1)=beta*intrate(+1)*xu6eas*(1+tceas)/(1+tceas(+1));
xu8eas(+1)*ggnam(+1)=beta*intrate(+1)*xu7eas*(1+tceas)/(1+tceas(+1));

zs2eas=phieas(-1)*((1-eduseas(-1))*(1-taueas(-1))*lams1eas(-1)*wagddseas(-1)+psieas(-1)*tr1seas(-1)*wagddseas(-1)+(1-eduseas(-1))*(1-lams1eas(-1))*bseas(-1)-(1+tceas(-1))*xs1eas(-1))/(ggnam*mmeas(-1));
zs3eas=intrate(-1)*zs2eas(-1)/(ggnam*mmeas(-1))+phieas(-2)*p2eas(-1)*((1-taueas(-1))*lams2eas(-1)*wagddseas(-1)+psieas(-1)*tr2seas(-1)*wagddseas(-1)+(1-lams2eas(-1))*bseas(-1)-(1+tceas(-1))*xs2eas(-1))/(ggnam*mmeas(-1)*mmeas(-2));
zs4eas=intrate(-1)*zs3eas(-1)/(ggnam*mmeas(-1))+phieas(-3)*p3eas(-1)*((1-taueas(-1))*lams3eas(-1)*wagddseas(-1)+psieas(-1)*tr3seas(-1)*wagddseas(-1)+(1-lams3eas(-1))*bseas(-1)-(1+tceas(-1))*xs3eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3));
zs5eas=intrate(-1)*zs4eas(-1)/(ggnam*mmeas(-1))+phieas(-4)*p4eas(-1)*((1-taueas(-1))*lams4eas(-1)*wagddseas(-1)+psieas(-1)*tr4seas(-1)*wagddseas(-1)+(1-lams4eas(-1))*bseas(-1)-(1+tceas(-1))*xs4eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4));
zs6eas=intrate(-1)*zs5eas(-1)/(ggnam*mmeas(-1))+phieas(-5)*p5eas(-1)*((1-taueas(-1))*lams5eas(-1)*wagddseas(-1)+psieas(-1)*tr5seas(-1)*wagddseas(-1)+(1-lams5eas(-1))*bseas(-1)-(1+tceas(-1))*xs5eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5));
zs7eas=intrate(-1)*zs6eas(-1)/(ggnam*mmeas(-1))+phieas(-6)*p6eas(-1)*((1-taueas(-1))*lams6eas(-1)*wagddseas(-1)+psieas(-1)*tr6seas(-1)*wagddseas(-1)+(1-lams6eas(-1))*bseas(-1)-(1+tceas(-1))*xs6eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6));
zs8eas=intrate(-1)*zs7eas(-1)/(ggnam*mmeas(-1))+phieas(-7)*p7eas(-1)*((1-taueas(-1))*lams7eas(-1)*wagddseas(-1)+psieas(-1)*tr7seas(-1)*wagddseas(-1)+(1-lams7eas(-1))*bseas(-1)-(1+tceas(-1))*xs7eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7));

zu2eas=(1-phieas(-1))*((1-eduueas(-1))*(1-taueas(-1))*lamu1eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr1ueas(-1)*wagddueas(-1)+(1-eduueas(-1))*(1-lamu1eas(-1))*bueas(-1)-(1+tceas(-1))*xu1eas(-1))/(ggnam*mmeas(-1));
zu3eas=intrate(-1)*zu2eas(-1)/(ggnam*mmeas(-1))+(1-phieas(-2))*p2eas(-1)*((1-taueas(-1))*lamu2eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr2ueas(-1)*wagddueas(-1)+(1-lamu2eas(-1))*bueas(-1)-(1+tceas(-1))*xu2eas(-1))/(ggnam*mmeas(-1)*mmeas(-2));
zu4eas=intrate(-1)*zu3eas(-1)/(ggnam*mmeas(-1))+(1-phieas(-3))*p3eas(-1)*((1-taueas(-1))*lamu3eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr3ueas(-1)*wagddueas(-1)+(1-lamu3eas(-1))*bueas(-1)-(1+tceas(-1))*xu3eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3));
zu5eas=intrate(-1)*zu4eas(-1)/(ggnam*mmeas(-1))+(1-phieas(-4))*p4eas(-1)*((1-taueas(-1))*lamu4eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr4ueas(-1)*wagddueas(-1)+(1-lamu4eas(-1))*bueas(-1)-(1+tceas(-1))*xu4eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4));
zu6eas=intrate(-1)*zu5eas(-1)/(ggnam*mmeas(-1))+(1-phieas(-5))*p5eas(-1)*((1-taueas(-1))*lamu5eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr5ueas(-1)*wagddueas(-1)+(1-lamu5eas(-1))*bueas(-1)-(1+tceas(-1))*xu5eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5));
zu7eas=intrate(-1)*zu6eas(-1)/(ggnam*mmeas(-1))+(1-phieas(-6))*p6eas(-1)*((1-taueas(-1))*lamu6eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr6ueas(-1)*wagddueas(-1)+(1-lamu6eas(-1))*bueas(-1)-(1+tceas(-1))*xu6eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6));
zu8eas=intrate(-1)*zu7eas(-1)/(ggnam*mmeas(-1))+(1-phieas(-7))*p7eas(-1)*((1-taueas(-1))*lamu7eas(-1)*(1-ureas(-1))*wagddueas(-1)+psieas(-1)*tr7ueas(-1)*wagddueas(-1)+(1-lamu7eas(-1))*bueas(-1)-(1+tceas(-1))*xu7eas(-1))/(ggnam*mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7));

weaeas=zs2eas+zs3eas+zs4eas+zs5eas+zs6eas+zs7eas+zs8eas
+zu2eas+zu3eas+zu4eas+zu5eas+zu6eas+zu7eas+zu8eas;

gdpeas=keas^alpha*(tfpeas*labeas)^(1-alpha);

goveas=taueas*(wagddseas*labseas+wagddueas*labddueas)
+pieas*intrate*keas
+tceas*(phieas*xs1eas
+xs2eas*phieas(-1)*p2eas/mmeas(-1)
+xs3eas*phieas(-2)*p3eas/(mmeas(-1)*mmeas(-2))
+xs4eas*phieas(-3)*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+xs5eas*phieas(-4)*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+xs6eas*phieas(-5)*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+xs7eas*phieas(-6)*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+xs8eas*phieas(-7)*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
+tceas*((1-phieas)*xu1eas
+xu2eas*(1-phieas(-1))*p2eas/mmeas(-1)
+xu3eas*(1-phieas(-2))*p3eas/(mmeas(-1)*mmeas(-2))
+xu4eas*(1-phieas(-3))*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+xu5eas*(1-phieas(-4))*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+xu6eas*(1-phieas(-5))*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+xu7eas*(1-phieas(-6))*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+xu8eas*(1-phieas(-7))*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
-bseas*((1-lams1eas)*phieas*(1-eduseas)
+(1-lams2eas)*phieas(-1)*p2eas/mmeas(-1)
+(1-lams3eas)*phieas(-2)*p3eas/(mmeas(-1)*mmeas(-2))
+(1-lams4eas)*phieas(-3)*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+(1-lams5eas)*phieas(-4)*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+(1-lams6eas)*phieas(-5)*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+(1-lams7eas)*phieas(-6)*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+(1-lams8eas)*phieas(-7)*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
-bueas*((1-lamu1eas)*(1-phieas)*(1-eduueas)
+(1-lamu2eas)*(1-phieas(-1))*p2eas/mmeas(-1)
+(1-lamu3eas)*(1-phieas(-2))*p3eas/(mmeas(-1)*mmeas(-2))
+(1-lamu4eas)*(1-phieas(-3))*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+(1-lamu5eas)*(1-phieas(-4))*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+(1-lamu6eas)*(1-phieas(-5))*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+(1-lamu7eas)*(1-phieas(-6))*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+(1-lamu8eas)*(1-phieas(-7))*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
-psieas*wagddseas*(tr1seas*phieas
+tr2seas*phieas(-1)*p2eas/mmeas(-1)
+tr3seas*phieas(-2)*p3eas/(mmeas(-1)*mmeas(-2))
+tr4seas*phieas(-3)*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+tr5seas*phieas(-4)*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+tr6seas*phieas(-5)*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+tr7seas*phieas(-6)*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+tr8seas*phieas(-7)*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
-psieas*wagddueas*(tr1ueas*(1-phieas)
+tr2ueas*(1-phieas(-1))*p2eas/mmeas(-1)
+tr3ueas*(1-phieas(-2))*p3eas/(mmeas(-1)*mmeas(-2))
+tr4ueas*(1-phieas(-3))*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+tr5ueas*(1-phieas(-4))*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+tr6ueas*(1-phieas(-5))*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+tr7ueas*(1-phieas(-6))*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+tr8ueas*(1-phieas(-7))*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
-conspubgdpeas
-pieas*intrate*keas;

// ------------------------------ china ------------------------------- //

xs8chi*(1+tcchi)=1/(p8chi*phichi(-7))*intrate*zs8chi*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)+lams8chi*(1-tauchi)*wagddschi+psichi*tr8schi*wagddschi+(1-lams8chi)*bschi;
xu8chi*(1+tcchi)=1/(p8chi*(1-phichi(-7)))*intrate*zu8chi*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)+(1-urchi)*lamu8chi*(1-tauchi)*wagdduchi+psichi*tr8uchi*wagdduchi+(1-lamu8chi)*buchi;

xs2chi(+1)*ggnam(+1)=beta*intrate(+1)*xs1chi*(1+tcchi)/(1+tcchi(+1));
xs3chi(+1)*ggnam(+1)=beta*intrate(+1)*xs2chi*(1+tcchi)/(1+tcchi(+1));
xs4chi(+1)*ggnam(+1)=beta*intrate(+1)*xs3chi*(1+tcchi)/(1+tcchi(+1));
xs5chi(+1)*ggnam(+1)=beta*intrate(+1)*xs4chi*(1+tcchi)/(1+tcchi(+1));
xs6chi(+1)*ggnam(+1)=beta*intrate(+1)*xs5chi*(1+tcchi)/(1+tcchi(+1));
xs7chi(+1)*ggnam(+1)=beta*intrate(+1)*xs6chi*(1+tcchi)/(1+tcchi(+1));
xs8chi(+1)*ggnam(+1)=beta*intrate(+1)*xs7chi*(1+tcchi)/(1+tcchi(+1));

xu2chi(+1)*ggnam(+1)=beta*intrate(+1)*xu1chi*(1+tcchi)/(1+tcchi(+1));
xu3chi(+1)*ggnam(+1)=beta*intrate(+1)*xu2chi*(1+tcchi)/(1+tcchi(+1));
xu4chi(+1)*ggnam(+1)=beta*intrate(+1)*xu3chi*(1+tcchi)/(1+tcchi(+1));
xu5chi(+1)*ggnam(+1)=beta*intrate(+1)*xu4chi*(1+tcchi)/(1+tcchi(+1));
xu6chi(+1)*ggnam(+1)=beta*intrate(+1)*xu5chi*(1+tcchi)/(1+tcchi(+1));
xu7chi(+1)*ggnam(+1)=beta*intrate(+1)*xu6chi*(1+tcchi)/(1+tcchi(+1));
xu8chi(+1)*ggnam(+1)=beta*intrate(+1)*xu7chi*(1+tcchi)/(1+tcchi(+1));

zs2chi=phichi(-1)*((1-eduschi(-1))*(1-tauchi(-1))*lams1chi(-1)*wagddschi(-1)+psichi(-1)*tr1schi(-1)*wagddschi(-1)+(1-eduschi(-1))*(1-lams1chi(-1))*bschi(-1)-(1+tcchi(-1))*xs1chi(-1))/(ggnam*mmchi(-1));
zs3chi=intrate(-1)*zs2chi(-1)/(ggnam*mmchi(-1))+phichi(-2)*p2chi(-1)*((1-tauchi(-1))*lams2chi(-1)*wagddschi(-1)+psichi(-1)*tr2schi(-1)*wagddschi(-1)+(1-lams2chi(-1))*bschi(-1)-(1+tcchi(-1))*xs2chi(-1))/(ggnam*mmchi(-1)*mmchi(-2));
zs4chi=intrate(-1)*zs3chi(-1)/(ggnam*mmchi(-1))+phichi(-3)*p3chi(-1)*((1-tauchi(-1))*lams3chi(-1)*wagddschi(-1)+psichi(-1)*tr3schi(-1)*wagddschi(-1)+(1-lams3chi(-1))*bschi(-1)-(1+tcchi(-1))*xs3chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3));
zs5chi=intrate(-1)*zs4chi(-1)/(ggnam*mmchi(-1))+phichi(-4)*p4chi(-1)*((1-tauchi(-1))*lams4chi(-1)*wagddschi(-1)+psichi(-1)*tr4schi(-1)*wagddschi(-1)+(1-lams4chi(-1))*bschi(-1)-(1+tcchi(-1))*xs4chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4));
zs6chi=intrate(-1)*zs5chi(-1)/(ggnam*mmchi(-1))+phichi(-5)*p5chi(-1)*((1-tauchi(-1))*lams5chi(-1)*wagddschi(-1)+psichi(-1)*tr5schi(-1)*wagddschi(-1)+(1-lams5chi(-1))*bschi(-1)-(1+tcchi(-1))*xs5chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5));
zs7chi=intrate(-1)*zs6chi(-1)/(ggnam*mmchi(-1))+phichi(-6)*p6chi(-1)*((1-tauchi(-1))*lams6chi(-1)*wagddschi(-1)+psichi(-1)*tr6schi(-1)*wagddschi(-1)+(1-lams6chi(-1))*bschi(-1)-(1+tcchi(-1))*xs6chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6));
zs8chi=intrate(-1)*zs7chi(-1)/(ggnam*mmchi(-1))+phichi(-7)*p7chi(-1)*((1-tauchi(-1))*lams7chi(-1)*wagddschi(-1)+psichi(-1)*tr7schi(-1)*wagddschi(-1)+(1-lams7chi(-1))*bschi(-1)-(1+tcchi(-1))*xs7chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7));

zu2chi=(1-phichi(-1))*((1-eduuchi(-1))*(1-tauchi(-1))*lamu1chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr1uchi(-1)*wagdduchi(-1)+(1-eduuchi(-1))*(1-lamu1chi(-1))*buchi(-1)-(1+tcchi(-1))*xu1chi(-1))/(ggnam*mmchi(-1));
zu3chi=intrate(-1)*zu2chi(-1)/(ggnam*mmchi(-1))+(1-phichi(-2))*p2chi(-1)*((1-tauchi(-1))*lamu2chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr2uchi(-1)*wagdduchi(-1)+(1-lamu2chi(-1))*buchi(-1)-(1+tcchi(-1))*xu2chi(-1))/(ggnam*mmchi(-1)*mmchi(-2));
zu4chi=intrate(-1)*zu3chi(-1)/(ggnam*mmchi(-1))+(1-phichi(-3))*p3chi(-1)*((1-tauchi(-1))*lamu3chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr3uchi(-1)*wagdduchi(-1)+(1-lamu3chi(-1))*buchi(-1)-(1+tcchi(-1))*xu3chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3));
zu5chi=intrate(-1)*zu4chi(-1)/(ggnam*mmchi(-1))+(1-phichi(-4))*p4chi(-1)*((1-tauchi(-1))*lamu4chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr4uchi(-1)*wagdduchi(-1)+(1-lamu4chi(-1))*buchi(-1)-(1+tcchi(-1))*xu4chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4));
zu6chi=intrate(-1)*zu5chi(-1)/(ggnam*mmchi(-1))+(1-phichi(-5))*p5chi(-1)*((1-tauchi(-1))*lamu5chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr5uchi(-1)*wagdduchi(-1)+(1-lamu5chi(-1))*buchi(-1)-(1+tcchi(-1))*xu5chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5));
zu7chi=intrate(-1)*zu6chi(-1)/(ggnam*mmchi(-1))+(1-phichi(-6))*p6chi(-1)*((1-tauchi(-1))*lamu6chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr6uchi(-1)*wagdduchi(-1)+(1-lamu6chi(-1))*buchi(-1)-(1+tcchi(-1))*xu6chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6));
zu8chi=intrate(-1)*zu7chi(-1)/(ggnam*mmchi(-1))+(1-phichi(-7))*p7chi(-1)*((1-tauchi(-1))*lamu7chi(-1)*(1-urchi(-1))*wagdduchi(-1)+psichi(-1)*tr7uchi(-1)*wagdduchi(-1)+(1-lamu7chi(-1))*buchi(-1)-(1+tcchi(-1))*xu7chi(-1))/(ggnam*mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7));

weachi=zs2chi+zs3chi+zs4chi+zs5chi+zs6chi+zs7chi+zs8chi
+zu2chi+zu3chi+zu4chi+zu5chi+zu6chi+zu7chi+zu8chi;

gdpchi=kchi^alpha*(tfpchi*labchi)^(1-alpha);

govchi=tauchi*(wagddschi*labschi+wagdduchi*labdduchi)
+pichi*intrate*kchi
+tcchi*(phichi*xs1chi
+xs2chi*phichi(-1)*p2chi/mmchi(-1)
+xs3chi*phichi(-2)*p3chi/(mmchi(-1)*mmchi(-2))
+xs4chi*phichi(-3)*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+xs5chi*phichi(-4)*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+xs6chi*phichi(-5)*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+xs7chi*phichi(-6)*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+xs8chi*phichi(-7)*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
+tcchi*((1-phichi)*xu1chi
+xu2chi*(1-phichi(-1))*p2chi/mmchi(-1)
+xu3chi*(1-phichi(-2))*p3chi/(mmchi(-1)*mmchi(-2))
+xu4chi*(1-phichi(-3))*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+xu5chi*(1-phichi(-4))*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+xu6chi*(1-phichi(-5))*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+xu7chi*(1-phichi(-6))*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+xu8chi*(1-phichi(-7))*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
-bschi*((1-lams1chi)*phichi*(1-eduschi)
+(1-lams2chi)*phichi(-1)*p2chi/mmchi(-1)
+(1-lams3chi)*phichi(-2)*p3chi/(mmchi(-1)*mmchi(-2))
+(1-lams4chi)*phichi(-3)*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+(1-lams5chi)*phichi(-4)*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+(1-lams6chi)*phichi(-5)*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+(1-lams7chi)*phichi(-6)*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+(1-lams8chi)*phichi(-7)*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
-buchi*((1-lamu1chi)*(1-phichi)*(1-eduuchi)
+(1-lamu2chi)*(1-phichi(-1))*p2chi/mmchi(-1)
+(1-lamu3chi)*(1-phichi(-2))*p3chi/(mmchi(-1)*mmchi(-2))
+(1-lamu4chi)*(1-phichi(-3))*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+(1-lamu5chi)*(1-phichi(-4))*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+(1-lamu6chi)*(1-phichi(-5))*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+(1-lamu7chi)*(1-phichi(-6))*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+(1-lamu8chi)*(1-phichi(-7))*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
-psichi*wagddschi*(tr1schi*phichi
+tr2schi*phichi(-1)*p2chi/mmchi(-1)
+tr3schi*phichi(-2)*p3chi/(mmchi(-1)*mmchi(-2))
+tr4schi*phichi(-3)*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+tr5schi*phichi(-4)*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+tr6schi*phichi(-5)*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+tr7schi*phichi(-6)*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+tr8schi*phichi(-7)*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
-psichi*wagdduchi*(tr1uchi*(1-phichi)
+tr2uchi*(1-phichi(-1))*p2chi/mmchi(-1)
+tr3uchi*(1-phichi(-2))*p3chi/(mmchi(-1)*mmchi(-2))
+tr4uchi*(1-phichi(-3))*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+tr5uchi*(1-phichi(-4))*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+tr6uchi*(1-phichi(-5))*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+tr7uchi*(1-phichi(-6))*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+tr8uchi*(1-phichi(-7))*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
-conspubgdpchi
-pichi*intrate*kchi;

// ------------------------------ india ------------------------------- //

xs8ind*(1+tcind)=1/(p8ind*phiind(-7))*intrate*zs8ind*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)+lams8ind*(1-tauind)*wagddsind+psiind*tr8sind*wagddsind+(1-lams8ind)*bsind;
xu8ind*(1+tcind)=1/(p8ind*(1-phiind(-7)))*intrate*zu8ind*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)+(1-urind)*lamu8ind*(1-tauind)*wagdduind+psiind*tr8uind*wagdduind+(1-lamu8ind)*buind;

xs2ind(+1)*ggnam(+1)=beta*intrate(+1)*xs1ind*(1+tcind)/(1+tcind(+1));
xs3ind(+1)*ggnam(+1)=beta*intrate(+1)*xs2ind*(1+tcind)/(1+tcind(+1));
xs4ind(+1)*ggnam(+1)=beta*intrate(+1)*xs3ind*(1+tcind)/(1+tcind(+1));
xs5ind(+1)*ggnam(+1)=beta*intrate(+1)*xs4ind*(1+tcind)/(1+tcind(+1));
xs6ind(+1)*ggnam(+1)=beta*intrate(+1)*xs5ind*(1+tcind)/(1+tcind(+1));
xs7ind(+1)*ggnam(+1)=beta*intrate(+1)*xs6ind*(1+tcind)/(1+tcind(+1));
xs8ind(+1)*ggnam(+1)=beta*intrate(+1)*xs7ind*(1+tcind)/(1+tcind(+1));

xu2ind(+1)*ggnam(+1)=beta*intrate(+1)*xu1ind*(1+tcind)/(1+tcind(+1));
xu3ind(+1)*ggnam(+1)=beta*intrate(+1)*xu2ind*(1+tcind)/(1+tcind(+1));
xu4ind(+1)*ggnam(+1)=beta*intrate(+1)*xu3ind*(1+tcind)/(1+tcind(+1));
xu5ind(+1)*ggnam(+1)=beta*intrate(+1)*xu4ind*(1+tcind)/(1+tcind(+1));
xu6ind(+1)*ggnam(+1)=beta*intrate(+1)*xu5ind*(1+tcind)/(1+tcind(+1));
xu7ind(+1)*ggnam(+1)=beta*intrate(+1)*xu6ind*(1+tcind)/(1+tcind(+1));
xu8ind(+1)*ggnam(+1)=beta*intrate(+1)*xu7ind*(1+tcind)/(1+tcind(+1));

zs2ind=phiind(-1)*((1-edusind(-1))*(1-tauind(-1))*lams1ind(-1)*wagddsind(-1)+psiind(-1)*tr1sind(-1)*wagddsind(-1)+(1-edusind(-1))*(1-lams1ind(-1))*bsind(-1)-(1+tcind(-1))*xs1ind(-1))/(ggnam*mmind(-1));
zs3ind=intrate(-1)*zs2ind(-1)/(ggnam*mmind(-1))+phiind(-2)*p2ind(-1)*((1-tauind(-1))*lams2ind(-1)*wagddsind(-1)+psiind(-1)*tr2sind(-1)*wagddsind(-1)+(1-lams2ind(-1))*bsind(-1)-(1+tcind(-1))*xs2ind(-1))/(ggnam*mmind(-1)*mmind(-2));
zs4ind=intrate(-1)*zs3ind(-1)/(ggnam*mmind(-1))+phiind(-3)*p3ind(-1)*((1-tauind(-1))*lams3ind(-1)*wagddsind(-1)+psiind(-1)*tr3sind(-1)*wagddsind(-1)+(1-lams3ind(-1))*bsind(-1)-(1+tcind(-1))*xs3ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3));
zs5ind=intrate(-1)*zs4ind(-1)/(ggnam*mmind(-1))+phiind(-4)*p4ind(-1)*((1-tauind(-1))*lams4ind(-1)*wagddsind(-1)+psiind(-1)*tr4sind(-1)*wagddsind(-1)+(1-lams4ind(-1))*bsind(-1)-(1+tcind(-1))*xs4ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4));
zs6ind=intrate(-1)*zs5ind(-1)/(ggnam*mmind(-1))+phiind(-5)*p5ind(-1)*((1-tauind(-1))*lams5ind(-1)*wagddsind(-1)+psiind(-1)*tr5sind(-1)*wagddsind(-1)+(1-lams5ind(-1))*bsind(-1)-(1+tcind(-1))*xs5ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5));
zs7ind=intrate(-1)*zs6ind(-1)/(ggnam*mmind(-1))+phiind(-6)*p6ind(-1)*((1-tauind(-1))*lams6ind(-1)*wagddsind(-1)+psiind(-1)*tr6sind(-1)*wagddsind(-1)+(1-lams6ind(-1))*bsind(-1)-(1+tcind(-1))*xs6ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6));
zs8ind=intrate(-1)*zs7ind(-1)/(ggnam*mmind(-1))+phiind(-7)*p7ind(-1)*((1-tauind(-1))*lams7ind(-1)*wagddsind(-1)+psiind(-1)*tr7sind(-1)*wagddsind(-1)+(1-lams7ind(-1))*bsind(-1)-(1+tcind(-1))*xs7ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7));

zu2ind=(1-phiind(-1))*((1-eduuind(-1))*(1-tauind(-1))*lamu1ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr1uind(-1)*wagdduind(-1)+(1-eduuind(-1))*(1-lamu1ind(-1))*buind(-1)-(1+tcind(-1))*xu1ind(-1))/(ggnam*mmind(-1));
zu3ind=intrate(-1)*zu2ind(-1)/(ggnam*mmind(-1))+(1-phiind(-2))*p2ind(-1)*((1-tauind(-1))*lamu2ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr2uind(-1)*wagdduind(-1)+(1-lamu2ind(-1))*buind(-1)-(1+tcind(-1))*xu2ind(-1))/(ggnam*mmind(-1)*mmind(-2));
zu4ind=intrate(-1)*zu3ind(-1)/(ggnam*mmind(-1))+(1-phiind(-3))*p3ind(-1)*((1-tauind(-1))*lamu3ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr3uind(-1)*wagdduind(-1)+(1-lamu3ind(-1))*buind(-1)-(1+tcind(-1))*xu3ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3));
zu5ind=intrate(-1)*zu4ind(-1)/(ggnam*mmind(-1))+(1-phiind(-4))*p4ind(-1)*((1-tauind(-1))*lamu4ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr4uind(-1)*wagdduind(-1)+(1-lamu4ind(-1))*buind(-1)-(1+tcind(-1))*xu4ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4));
zu6ind=intrate(-1)*zu5ind(-1)/(ggnam*mmind(-1))+(1-phiind(-5))*p5ind(-1)*((1-tauind(-1))*lamu5ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr5uind(-1)*wagdduind(-1)+(1-lamu5ind(-1))*buind(-1)-(1+tcind(-1))*xu5ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5));
zu7ind=intrate(-1)*zu6ind(-1)/(ggnam*mmind(-1))+(1-phiind(-6))*p6ind(-1)*((1-tauind(-1))*lamu6ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr6uind(-1)*wagdduind(-1)+(1-lamu6ind(-1))*buind(-1)-(1+tcind(-1))*xu6ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6));
zu8ind=intrate(-1)*zu7ind(-1)/(ggnam*mmind(-1))+(1-phiind(-7))*p7ind(-1)*((1-tauind(-1))*lamu7ind(-1)*(1-urind(-1))*wagdduind(-1)+psiind(-1)*tr7uind(-1)*wagdduind(-1)+(1-lamu7ind(-1))*buind(-1)-(1+tcind(-1))*xu7ind(-1))/(ggnam*mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7));

weaind=zs2ind+zs3ind+zs4ind+zs5ind+zs6ind+zs7ind+zs8ind
+zu2ind+zu3ind+zu4ind+zu5ind+zu6ind+zu7ind+zu8ind;

gdpind=kind^alpha*(tfpind*labind)^(1-alpha);

govind=tauind*(wagddsind*labsind+wagdduind*labdduind)
+piind*intrate*kind
+tcind*(phiind*xs1ind
+xs2ind*phiind(-1)*p2ind/mmind(-1)
+xs3ind*phiind(-2)*p3ind/(mmind(-1)*mmind(-2))
+xs4ind*phiind(-3)*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+xs5ind*phiind(-4)*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+xs6ind*phiind(-5)*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+xs7ind*phiind(-6)*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+xs8ind*phiind(-7)*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
+tcind*((1-phiind)*xu1ind
+xu2ind*(1-phiind(-1))*p2ind/mmind(-1)
+xu3ind*(1-phiind(-2))*p3ind/(mmind(-1)*mmind(-2))
+xu4ind*(1-phiind(-3))*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+xu5ind*(1-phiind(-4))*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+xu6ind*(1-phiind(-5))*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+xu7ind*(1-phiind(-6))*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+xu8ind*(1-phiind(-7))*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
-bsind*((1-lams1ind)*phiind*(1-edusind)
+(1-lams2ind)*phiind(-1)*p2ind/mmind(-1)
+(1-lams3ind)*phiind(-2)*p3ind/(mmind(-1)*mmind(-2))
+(1-lams4ind)*phiind(-3)*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+(1-lams5ind)*phiind(-4)*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+(1-lams6ind)*phiind(-5)*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+(1-lams7ind)*phiind(-6)*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+(1-lams8ind)*phiind(-7)*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
-buind*((1-lamu1ind)*(1-phiind)*(1-eduuind)
+(1-lamu2ind)*(1-phiind(-1))*p2ind/mmind(-1)
+(1-lamu3ind)*(1-phiind(-2))*p3ind/(mmind(-1)*mmind(-2))
+(1-lamu4ind)*(1-phiind(-3))*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+(1-lamu5ind)*(1-phiind(-4))*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+(1-lamu6ind)*(1-phiind(-5))*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+(1-lamu7ind)*(1-phiind(-6))*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+(1-lamu8ind)*(1-phiind(-7))*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
-psiind*wagddsind*(tr1sind*phiind
+tr2sind*phiind(-1)*p2ind/mmind(-1)
+tr3sind*phiind(-2)*p3ind/(mmind(-1)*mmind(-2))
+tr4sind*phiind(-3)*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+tr5sind*phiind(-4)*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+tr6sind*phiind(-5)*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+tr7sind*phiind(-6)*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+tr8sind*phiind(-7)*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
-psiind*wagdduind*(tr1uind*(1-phiind)
+tr2uind*(1-phiind(-1))*p2ind/mmind(-1)
+tr3uind*(1-phiind(-2))*p3ind/(mmind(-1)*mmind(-2))
+tr4uind*(1-phiind(-3))*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+tr5uind*(1-phiind(-4))*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+tr6uind*(1-phiind(-5))*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+tr7uind*(1-phiind(-6))*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+tr8uind*(1-phiind(-7))*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
-conspubgdpind
-piind*intrate*kind;

//--------------------------- bengdp --------------------------------- //

bengdpadv=(bsadv*((1-lams1adv)*phiadv*(1-edusadv)
+(1-lams2adv)*phiadv(-1)*p2adv/mmadv(-1)
+(1-lams3adv)*phiadv(-2)*p3adv/(mmadv(-1)*mmadv(-2))
+(1-lams4adv)*phiadv(-3)*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+(1-lams5adv)*phiadv(-4)*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+(1-lams6adv)*phiadv(-5)*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+(1-lams7adv)*phiadv(-6)*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+(1-lams8adv)*phiadv(-7)*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
+buadv*((1-lamu1adv)*(1-phiadv)*(1-eduuadv)
+(1-lamu2adv)*(1-phiadv(-1))*p2adv/mmadv(-1)
+(1-lamu3adv)*(1-phiadv(-2))*p3adv/(mmadv(-1)*mmadv(-2))
+(1-lamu4adv)*(1-phiadv(-3))*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+(1-lamu5adv)*(1-phiadv(-4))*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+(1-lamu6adv)*(1-phiadv(-5))*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+(1-lamu7adv)*(1-phiadv(-6))*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+(1-lamu8adv)*(1-phiadv(-7))*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7))))/gdpadv;

bengdpnam=(bsnam*((1-lams1nam)*phinam*(1-edusnam)
+(1-lams2nam)*phinam(-1)*p2nam/mmnam(-1)
+(1-lams3nam)*phinam(-2)*p3nam/(mmnam(-1)*mmnam(-2))
+(1-lams4nam)*phinam(-3)*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+(1-lams5nam)*phinam(-4)*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+(1-lams6nam)*phinam(-5)*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+(1-lams7nam)*phinam(-6)*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+(1-lams8nam)*phinam(-7)*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
+bunam*((1-lamu1nam)*(1-phinam)*(1-eduunam)
+(1-lamu2nam)*(1-phinam(-1))*p2nam/mmnam(-1)
+(1-lamu3nam)*(1-phinam(-2))*p3nam/(mmnam(-1)*mmnam(-2))
+(1-lamu4nam)*(1-phinam(-3))*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+(1-lamu5nam)*(1-phinam(-4))*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+(1-lamu6nam)*(1-phinam(-5))*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+(1-lamu7nam)*(1-phinam(-6))*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+(1-lamu8nam)*(1-phinam(-7))*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7))))/gdpnam;

bengdpssa=(bsssa*((1-lams1ssa)*phissa*(1-edusssa)
+(1-lams2ssa)*phissa(-1)*p2ssa/mmssa(-1)
+(1-lams3ssa)*phissa(-2)*p3ssa/(mmssa(-1)*mmssa(-2))
+(1-lams4ssa)*phissa(-3)*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+(1-lams5ssa)*phissa(-4)*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+(1-lams6ssa)*phissa(-5)*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+(1-lams7ssa)*phissa(-6)*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+(1-lams8ssa)*phissa(-7)*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
+bussa*((1-lamu1ssa)*(1-phissa)*(1-eduussa)
+(1-lamu2ssa)*(1-phissa(-1))*p2ssa/mmssa(-1)
+(1-lamu3ssa)*(1-phissa(-2))*p3ssa/(mmssa(-1)*mmssa(-2))
+(1-lamu4ssa)*(1-phissa(-3))*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+(1-lamu5ssa)*(1-phissa(-4))*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+(1-lamu6ssa)*(1-phissa(-5))*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+(1-lamu7ssa)*(1-phissa(-6))*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+(1-lamu8ssa)*(1-phissa(-7))*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7))))/gdpssa;

bengdplac=(bslac*((1-lams1lac)*philac*(1-eduslac)
+(1-lams2lac)*philac(-1)*p2lac/mmlac(-1)
+(1-lams3lac)*philac(-2)*p3lac/(mmlac(-1)*mmlac(-2))
+(1-lams4lac)*philac(-3)*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+(1-lams5lac)*philac(-4)*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+(1-lams6lac)*philac(-5)*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+(1-lams7lac)*philac(-6)*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+(1-lams8lac)*philac(-7)*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
+bulac*((1-lamu1lac)*(1-philac)*(1-eduulac)
+(1-lamu2lac)*(1-philac(-1))*p2lac/mmlac(-1)
+(1-lamu3lac)*(1-philac(-2))*p3lac/(mmlac(-1)*mmlac(-2))
+(1-lamu4lac)*(1-philac(-3))*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+(1-lamu5lac)*(1-philac(-4))*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+(1-lamu6lac)*(1-philac(-5))*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+(1-lamu7lac)*(1-philac(-6))*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+(1-lamu8lac)*(1-philac(-7))*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7))))/gdplac;

bengdpjap=(bsjap*((1-lams1jap)*phijap*(1-edusjap)
+(1-lams2jap)*phijap(-1)*p2jap/mmjap(-1)
+(1-lams3jap)*phijap(-2)*p3jap/(mmjap(-1)*mmjap(-2))
+(1-lams4jap)*phijap(-3)*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+(1-lams5jap)*phijap(-4)*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+(1-lams6jap)*phijap(-5)*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+(1-lams7jap)*phijap(-6)*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+(1-lams8jap)*phijap(-7)*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
+bujap*((1-lamu1jap)*(1-phijap)*(1-eduujap)
+(1-lamu2jap)*(1-phijap(-1))*p2jap/mmjap(-1)
+(1-lamu3jap)*(1-phijap(-2))*p3jap/(mmjap(-1)*mmjap(-2))
+(1-lamu4jap)*(1-phijap(-3))*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+(1-lamu5jap)*(1-phijap(-4))*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+(1-lamu6jap)*(1-phijap(-5))*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+(1-lamu7jap)*(1-phijap(-6))*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+(1-lamu8jap)*(1-phijap(-7))*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7))))/gdpjap;

bengdpeas=(bseas*((1-lams1eas)*phieas*(1-eduseas)
+(1-lams2eas)*phieas(-1)*p2eas/mmeas(-1)
+(1-lams3eas)*phieas(-2)*p3eas/(mmeas(-1)*mmeas(-2))
+(1-lams4eas)*phieas(-3)*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+(1-lams5eas)*phieas(-4)*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+(1-lams6eas)*phieas(-5)*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+(1-lams7eas)*phieas(-6)*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+(1-lams8eas)*phieas(-7)*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
+bueas*((1-lamu1eas)*(1-phieas)*(1-eduueas)
+(1-lamu2eas)*(1-phieas(-1))*p2eas/mmeas(-1)
+(1-lamu3eas)*(1-phieas(-2))*p3eas/(mmeas(-1)*mmeas(-2))
+(1-lamu4eas)*(1-phieas(-3))*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+(1-lamu5eas)*(1-phieas(-4))*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+(1-lamu6eas)*(1-phieas(-5))*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+(1-lamu7eas)*(1-phieas(-6))*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+(1-lamu8eas)*(1-phieas(-7))*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7))))/gdpeas;

bengdprus=(bsrus*((1-lams1rus)*phirus*(1-edusrus)
+(1-lams2rus)*phirus(-1)*p2rus/mmrus(-1)
+(1-lams3rus)*phirus(-2)*p3rus/(mmrus(-1)*mmrus(-2))
+(1-lams4rus)*phirus(-3)*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+(1-lams5rus)*phirus(-4)*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+(1-lams6rus)*phirus(-5)*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+(1-lams7rus)*phirus(-6)*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+(1-lams8rus)*phirus(-7)*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
+burus*((1-lamu1rus)*(1-phirus)*(1-eduurus)
+(1-lamu2rus)*(1-phirus(-1))*p2rus/mmrus(-1)
+(1-lamu3rus)*(1-phirus(-2))*p3rus/(mmrus(-1)*mmrus(-2))
+(1-lamu4rus)*(1-phirus(-3))*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+(1-lamu5rus)*(1-phirus(-4))*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+(1-lamu6rus)*(1-phirus(-5))*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+(1-lamu7rus)*(1-phirus(-6))*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+(1-lamu8rus)*(1-phirus(-7))*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7))))/gdprus;

bengdpmen=(bsmen*((1-lams1men)*phimen*(1-edusmen)
+(1-lams2men)*phimen(-1)*p2men/mmmen(-1)
+(1-lams3men)*phimen(-2)*p3men/(mmmen(-1)*mmmen(-2))
+(1-lams4men)*phimen(-3)*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+(1-lams5men)*phimen(-4)*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+(1-lams6men)*phimen(-5)*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+(1-lams7men)*phimen(-6)*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+(1-lams8men)*phimen(-7)*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
+bumen*((1-lamu1men)*(1-phimen)*(1-eduumen)
+(1-lamu2men)*(1-phimen(-1))*p2men/mmmen(-1)
+(1-lamu3men)*(1-phimen(-2))*p3men/(mmmen(-1)*mmmen(-2))
+(1-lamu4men)*(1-phimen(-3))*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+(1-lamu5men)*(1-phimen(-4))*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+(1-lamu6men)*(1-phimen(-5))*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+(1-lamu7men)*(1-phimen(-6))*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+(1-lamu8men)*(1-phimen(-7))*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7))))/gdpmen;

bengdpchi=(bschi*((1-lams1chi)*phichi*(1-eduschi)
+(1-lams2chi)*phichi(-1)*p2chi/mmchi(-1)
+(1-lams3chi)*phichi(-2)*p3chi/(mmchi(-1)*mmchi(-2))
+(1-lams4chi)*phichi(-3)*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+(1-lams5chi)*phichi(-4)*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+(1-lams6chi)*phichi(-5)*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+(1-lams7chi)*phichi(-6)*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+(1-lams8chi)*phichi(-7)*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
+buchi*((1-lamu1chi)*(1-phichi)*(1-eduuchi)
+(1-lamu2chi)*(1-phichi(-1))*p2chi/mmchi(-1)
+(1-lamu3chi)*(1-phichi(-2))*p3chi/(mmchi(-1)*mmchi(-2))
+(1-lamu4chi)*(1-phichi(-3))*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+(1-lamu5chi)*(1-phichi(-4))*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+(1-lamu6chi)*(1-phichi(-5))*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+(1-lamu7chi)*(1-phichi(-6))*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+(1-lamu8chi)*(1-phichi(-7))*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7))))/gdpchi;

bengdpind=(bsind*((1-lams1ind)*phiind*(1-edusind)
+(1-lams2ind)*phiind(-1)*p2ind/mmind(-1)
+(1-lams3ind)*phiind(-2)*p3ind/(mmind(-1)*mmind(-2))
+(1-lams4ind)*phiind(-3)*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+(1-lams5ind)*phiind(-4)*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+(1-lams6ind)*phiind(-5)*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+(1-lams7ind)*phiind(-6)*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+(1-lams8ind)*phiind(-7)*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
+buind*((1-lamu1ind)*(1-phiind)*(1-eduuind)
+(1-lamu2ind)*(1-phiind(-1))*p2ind/mmind(-1)
+(1-lamu3ind)*(1-phiind(-2))*p3ind/(mmind(-1)*mmind(-2))
+(1-lamu4ind)*(1-phiind(-3))*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+(1-lamu5ind)*(1-phiind(-4))*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+(1-lamu6ind)*(1-phiind(-5))*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+(1-lamu7ind)*(1-phiind(-6))*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+(1-lamu8ind)*(1-phiind(-7))*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7))))/gdpind;

// ----------------------------- sup ----------------------------------- //

supadv=(1+p2adv/mmadv(-1)
+p3adv/(mmadv(-1)*mmadv(-2))
+p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)))/
(p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)));

supnam=(1+p2nam/mmnam(-1)
+p3nam/(mmnam(-1)*mmnam(-2))
+p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)))/
(p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

supssa=(1+p2ssa/mmssa(-1)
+p3ssa/(mmssa(-1)*mmssa(-2))
+p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)))/
(p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)));

suplac=(1+p2lac/mmlac(-1)
+p3lac/(mmlac(-1)*mmlac(-2))
+p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)))/
(p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)));

supeas=(1+p2eas/mmeas(-1)
+p3eas/(mmeas(-1)*mmeas(-2))
+p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)))/
(p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)));

supmen=(1+p2men/mmmen(-1)
+p3men/(mmmen(-1)*mmmen(-2))
+p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)))/
(p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)));

suprus=(1+p2rus/mmrus(-1)
+p3rus/(mmrus(-1)*mmrus(-2))
+p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)))/
(p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)));

supjap=(1+p2jap/mmjap(-1)
+p3jap/(mmjap(-1)*mmjap(-2))
+p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)))/
(p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)));

supchi=(1+p2chi/mmchi(-1)
+p3chi/(mmchi(-1)*mmchi(-2))
+p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)))/
(p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)));

supind=(1+p2ind/mmind(-1)
+p3ind/(mmind(-1)*mmind(-2))
+p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)))/
(p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)));

// ----------------------------- trgdp -------------------------------- //

trgdpadv=(psiadv*wagddsadv*(tr1sadv*phiadv
+tr2sadv*phiadv(-1)*p2adv/mmadv(-1)
+tr3sadv*phiadv(-2)*p3adv/(mmadv(-1)*mmadv(-2))
+tr4sadv*phiadv(-3)*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+tr5sadv*phiadv(-4)*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+tr6sadv*phiadv(-5)*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+tr7sadv*phiadv(-6)*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+tr8sadv*phiadv(-7)*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
+psiadv*wagdduadv*(tr1uadv*(1-phiadv)
+tr2uadv*(1-phiadv(-1))*p2adv/mmadv(-1)
+tr3uadv*(1-phiadv(-2))*p3adv/(mmadv(-1)*mmadv(-2))
+tr4uadv*(1-phiadv(-3))*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+tr5uadv*(1-phiadv(-4))*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+tr6uadv*(1-phiadv(-5))*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+tr7uadv*(1-phiadv(-6))*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+tr8uadv*(1-phiadv(-7))*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7))))/gdpadv;

trgdpnam=(psinam*wagddsnam*(tr1snam*phinam
+tr2snam*phinam(-1)*p2nam/mmnam(-1)
+tr3snam*phinam(-2)*p3nam/(mmnam(-1)*mmnam(-2))
+tr4snam*phinam(-3)*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+tr5snam*phinam(-4)*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+tr6snam*phinam(-5)*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+tr7snam*phinam(-6)*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+tr8snam*phinam(-7)*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
+psinam*wagddunam*(tr1unam*(1-phinam)
+tr2unam*(1-phinam(-1))*p2nam/mmnam(-1)
+tr3unam*(1-phinam(-2))*p3nam/(mmnam(-1)*mmnam(-2))
+tr4unam*(1-phinam(-3))*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+tr5unam*(1-phinam(-4))*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+tr6unam*(1-phinam(-5))*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+tr7unam*(1-phinam(-6))*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+tr8unam*(1-phinam(-7))*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7))))/gdpnam;

trgdpssa=(psissa*wagddsssa*(tr1sssa*phissa
+tr2sssa*phissa(-1)*p2ssa/mmssa(-1)
+tr3sssa*phissa(-2)*p3ssa/(mmssa(-1)*mmssa(-2))
+tr4sssa*phissa(-3)*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+tr5sssa*phissa(-4)*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+tr6sssa*phissa(-5)*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+tr7sssa*phissa(-6)*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+tr8sssa*phissa(-7)*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
+psissa*wagddussa*(tr1ussa*(1-phissa)
+tr2ussa*(1-phissa(-1))*p2ssa/mmssa(-1)
+tr3ussa*(1-phissa(-2))*p3ssa/(mmssa(-1)*mmssa(-2))
+tr4ussa*(1-phissa(-3))*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+tr5ussa*(1-phissa(-4))*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+tr6ussa*(1-phissa(-5))*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+tr7ussa*(1-phissa(-6))*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+tr8ussa*(1-phissa(-7))*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7))))/gdpssa;

trgdplac=(psilac*wagddslac*(tr1slac*philac
+tr2slac*philac(-1)*p2lac/mmlac(-1)
+tr3slac*philac(-2)*p3lac/(mmlac(-1)*mmlac(-2))
+tr4slac*philac(-3)*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+tr5slac*philac(-4)*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+tr6slac*philac(-5)*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+tr7slac*philac(-6)*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+tr8slac*philac(-7)*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
+psilac*wagddulac*(tr1ulac*(1-philac)
+tr2ulac*(1-philac(-1))*p2lac/mmlac(-1)
+tr3ulac*(1-philac(-2))*p3lac/(mmlac(-1)*mmlac(-2))
+tr4ulac*(1-philac(-3))*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+tr5ulac*(1-philac(-4))*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+tr6ulac*(1-philac(-5))*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+tr7ulac*(1-philac(-6))*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+tr8ulac*(1-philac(-7))*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7))))/gdplac;

trgdpjap=(psijap*hjap*wagddsjap*(tr1sjap*phijap
+tr2sjap*phijap(-1)*p2jap/mmjap(-1)
+tr3sjap*phijap(-2)*p3jap/(mmjap(-1)*mmjap(-2))
+tr4sjap*phijap(-3)*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+tr5sjap*phijap(-4)*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+tr6sjap*phijap(-5)*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+tr7sjap*phijap(-6)*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+tr8sjap*phijap(-7)*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
+psijap*wagddujap*(tr1ujap*(1-phijap)
+tr2ujap*(1-phijap(-1))*p2jap/mmjap(-1)
+tr3ujap*(1-phijap(-2))*p3jap/(mmjap(-1)*mmjap(-2))
+tr4ujap*(1-phijap(-3))*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+tr5ujap*(1-phijap(-4))*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+tr6ujap*(1-phijap(-5))*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+tr7ujap*(1-phijap(-6))*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+tr8ujap*(1-phijap(-7))*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7))))/gdpjap;

trgdprus=(psirus*wagddsrus*(tr1srus*phirus
+tr2srus*phirus(-1)*p2rus/mmrus(-1)
+tr3srus*phirus(-2)*p3rus/(mmrus(-1)*mmrus(-2))
+tr4srus*phirus(-3)*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+tr5srus*phirus(-4)*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+tr6srus*phirus(-5)*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+tr7srus*phirus(-6)*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+tr8srus*phirus(-7)*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
+psirus*wagddurus*(tr1urus*(1-phirus)
+tr2urus*(1-phirus(-1))*p2rus/mmrus(-1)
+tr3urus*(1-phirus(-2))*p3rus/(mmrus(-1)*mmrus(-2))
+tr4urus*(1-phirus(-3))*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+tr5urus*(1-phirus(-4))*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+tr6urus*(1-phirus(-5))*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+tr7urus*(1-phirus(-6))*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+tr8urus*(1-phirus(-7))*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7))))/gdprus;

trgdpmen=(psimen*wagddsmen*(tr1smen*phimen
+tr2smen*phimen(-1)*p2men/mmmen(-1)
+tr3smen*phimen(-2)*p3men/(mmmen(-1)*mmmen(-2))
+tr4smen*phimen(-3)*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+tr5smen*phimen(-4)*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+tr6smen*phimen(-5)*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+tr7smen*phimen(-6)*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+tr8smen*phimen(-7)*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
+psimen*wagddumen*(tr1umen*(1-phimen)
+tr2umen*(1-phimen(-1))*p2men/mmmen(-1)
+tr3umen*(1-phimen(-2))*p3men/(mmmen(-1)*mmmen(-2))
+tr4umen*(1-phimen(-3))*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+tr5umen*(1-phimen(-4))*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+tr6umen*(1-phimen(-5))*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+tr7umen*(1-phimen(-6))*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+tr8umen*(1-phimen(-7))*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7))))/gdpmen;

trgdpeas=(psieas*wagddseas*(tr1seas*phieas
+tr2seas*phieas(-1)*p2eas/mmeas(-1)
+tr3seas*phieas(-2)*p3eas/(mmeas(-1)*mmeas(-2))
+tr4seas*phieas(-3)*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+tr5seas*phieas(-4)*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+tr6seas*phieas(-5)*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+tr7seas*phieas(-6)*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+tr8seas*phieas(-7)*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
+psieas*wagddueas*(tr1ueas*(1-phieas)
+tr2ueas*(1-phieas(-1))*p2eas/mmeas(-1)
+tr3ueas*(1-phieas(-2))*p3eas/(mmeas(-1)*mmeas(-2))
+tr4ueas*(1-phieas(-3))*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+tr5ueas*(1-phieas(-4))*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+tr6ueas*(1-phieas(-5))*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+tr7ueas*(1-phieas(-6))*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+tr8ueas*(1-phieas(-7))*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7))))/gdpeas;

trgdpchi=(psichi*wagddschi*(tr1schi*phichi
+tr2schi*phichi(-1)*p2chi/mmchi(-1)
+tr3schi*phichi(-2)*p3chi/(mmchi(-1)*mmchi(-2))
+tr4schi*phichi(-3)*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+tr5schi*phichi(-4)*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+tr6schi*phichi(-5)*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+tr7schi*phichi(-6)*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+tr8schi*phichi(-7)*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
+psichi*wagdduchi*(tr1uchi*(1-phichi)
+tr2uchi*(1-phichi(-1))*p2chi/mmchi(-1)
+tr3uchi*(1-phichi(-2))*p3chi/(mmchi(-1)*mmchi(-2))
+tr4uchi*(1-phichi(-3))*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+tr5uchi*(1-phichi(-4))*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+tr6uchi*(1-phichi(-5))*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+tr7uchi*(1-phichi(-6))*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+tr8uchi*(1-phichi(-7))*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7))))/gdpchi;

trgdpind=(psiind*wagddsind*(tr1sind*phiind
+tr2sind*phiind(-1)*p2ind/mmind(-1)
+tr3sind*phiind(-2)*p3ind/(mmind(-1)*mmind(-2))
+tr4sind*phiind(-3)*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+tr5sind*phiind(-4)*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+tr6sind*phiind(-5)*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+tr7sind*phiind(-6)*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+tr8sind*phiind(-7)*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
+psiind*wagdduind*(tr1uind*(1-phiind)
+tr2uind*(1-phiind(-1))*p2ind/mmind(-1)
+tr3uind*(1-phiind(-2))*p3ind/(mmind(-1)*mmind(-2))
+tr4uind*(1-phiind(-3))*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+tr5uind*(1-phiind(-4))*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+tr6uind*(1-phiind(-5))*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+tr7uind*(1-phiind(-6))*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+tr8uind*(1-phiind(-7))*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7))))/gdpind;

// ------------------------- rapport lab-lab sans h ------------------- //

lablabsshadv=labadv/((1-edusadv)*phiadv*lams1adv
+phiadv(-1)*lams2adv*p2adv/mmadv(-1)
+phiadv(-2)*lams3adv*p3adv/(mmadv(-1)*mmadv(-2))
+phiadv(-3)*lams4adv*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+phiadv(-4)*lams5adv*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+phiadv(-5)*lams6adv*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+phiadv(-6)*lams7adv*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+phiadv(-7)*lams8adv*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7))
+(1-eduuadv)*(1-phiadv)*lamu1adv
+(1-phiadv(-1))*lamu2adv*p2adv/mmadv(-1)
+(1-phiadv(-2))*lamu3adv*p3adv/(mmadv(-1)*mmadv(-2))
+(1-phiadv(-3))*lamu4adv*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+(1-phiadv(-4))*lamu5adv*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+(1-phiadv(-5))*lamu6adv*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+(1-phiadv(-6))*lamu7adv*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+(1-phiadv(-7))*lamu8adv*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)));

lablabsshnam=labnam/((1-edusnam)*phinam*lams1nam
+phinam(-1)*lams2nam*p2nam/mmnam(-1)
+phinam(-2)*lams3nam*p3nam/(mmnam(-1)*mmnam(-2))
+phinam(-3)*lams4nam*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+phinam(-4)*lams5nam*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+phinam(-5)*lams6nam*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+phinam(-6)*lams7nam*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+phinam(-7)*lams8nam*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7))
+(1-eduunam)*(1-phinam)*lamu1nam
+(1-phinam(-1))*lamu2nam*p2nam/mmnam(-1)
+(1-phinam(-2))*lamu3nam*p3nam/(mmnam(-1)*mmnam(-2))
+(1-phinam(-3))*lamu4nam*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+(1-phinam(-4))*lamu5nam*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+(1-phinam(-5))*lamu6nam*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+(1-phinam(-6))*lamu7nam*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+(1-phinam(-7))*lamu8nam*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)));

lablabsshjap=labjap/((1-edusjap)*phijap*lams1jap
+phijap(-1)*lams2jap*p2jap/mmjap(-1)
+phijap(-2)*lams3jap*p3jap/(mmjap(-1)*mmjap(-2))
+phijap(-3)*lams4jap*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+phijap(-4)*lams5jap*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+phijap(-5)*lams6jap*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+phijap(-6)*lams7jap*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+phijap(-7)*lams8jap*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7))
+(1-eduujap)*(1-phijap)*lamu1jap
+(1-phijap(-1))*lamu2jap*p2jap/mmjap(-1)
+(1-phijap(-2))*lamu3jap*p3jap/(mmjap(-1)*mmjap(-2))
+(1-phijap(-3))*lamu4jap*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+(1-phijap(-4))*lamu5jap*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+(1-phijap(-5))*lamu6jap*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+(1-phijap(-6))*lamu7jap*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+(1-phijap(-7))*lamu8jap*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)));

lablabsshlac=lablac/((1-eduslac)*philac*lams1lac
+philac(-1)*lams2lac*p2lac/mmlac(-1)
+philac(-2)*lams3lac*p3lac/(mmlac(-1)*mmlac(-2))
+philac(-3)*lams4lac*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+philac(-4)*lams5lac*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+philac(-5)*lams6lac*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+philac(-6)*lams7lac*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+philac(-7)*lams8lac*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7))
+(1-eduulac)*(1-philac)*lamu1lac
+(1-philac(-1))*lamu2lac*p2lac/mmlac(-1)
+(1-philac(-2))*lamu3lac*p3lac/(mmlac(-1)*mmlac(-2))
+(1-philac(-3))*lamu4lac*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+(1-philac(-4))*lamu5lac*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+(1-philac(-5))*lamu6lac*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+(1-philac(-6))*lamu7lac*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+(1-philac(-7))*lamu8lac*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)));

lablabsshssa=labssa/((1-edusssa)*phissa*lams1ssa
+phissa(-1)*lams2ssa*p2ssa/mmssa(-1)
+phissa(-2)*lams3ssa*p3ssa/(mmssa(-1)*mmssa(-2))
+phissa(-3)*lams4ssa*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+phissa(-4)*lams5ssa*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+phissa(-5)*lams6ssa*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+phissa(-6)*lams7ssa*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+phissa(-7)*lams8ssa*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7))
+(1-eduussa)*(1-phissa)*lamu1ssa
+(1-phissa(-1))*lamu2ssa*p2ssa/mmssa(-1)
+(1-phissa(-2))*lamu3ssa*p3ssa/(mmssa(-1)*mmssa(-2))
+(1-phissa(-3))*lamu4ssa*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+(1-phissa(-4))*lamu5ssa*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+(1-phissa(-5))*lamu6ssa*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+(1-phissa(-6))*lamu7ssa*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+(1-phissa(-7))*lamu8ssa*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)));

lablabssheas=labeas/((1-eduseas)*phieas*lams1eas
+phieas(-1)*lams2eas*p2eas/mmeas(-1)
+phieas(-2)*lams3eas*p3eas/(mmeas(-1)*mmeas(-2))
+phieas(-3)*lams4eas*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+phieas(-4)*lams5eas*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+phieas(-5)*lams6eas*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+phieas(-6)*lams7eas*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+phieas(-7)*lams8eas*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7))
+(1-eduueas)*(1-phieas)*lamu1eas
+(1-phieas(-1))*lamu2eas*p2eas/mmeas(-1)
+(1-phieas(-2))*lamu3eas*p3eas/(mmeas(-1)*mmeas(-2))
+(1-phieas(-3))*lamu4eas*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+(1-phieas(-4))*lamu5eas*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+(1-phieas(-5))*lamu6eas*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+(1-phieas(-6))*lamu7eas*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+(1-phieas(-7))*lamu8eas*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)));

lablabsshmen=labmen/((1-edusmen)*phimen*lams1men
+phimen(-1)*lams2men*p2men/mmmen(-1)
+phimen(-2)*lams3men*p3men/(mmmen(-1)*mmmen(-2))
+phimen(-3)*lams4men*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+phimen(-4)*lams5men*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+phimen(-5)*lams6men*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+phimen(-6)*lams7men*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+phimen(-7)*lams8men*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7))
+(1-eduumen)*(1-phimen)*lamu1men
+(1-phimen(-1))*lamu2men*p2men/mmmen(-1)
+(1-phimen(-2))*lamu3men*p3men/(mmmen(-1)*mmmen(-2))
+(1-phimen(-3))*lamu4men*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+(1-phimen(-4))*lamu5men*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+(1-phimen(-5))*lamu6men*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+(1-phimen(-6))*lamu7men*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+(1-phimen(-7))*lamu8men*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)));

lablabsshrus=labrus/((1-edusrus)*phirus*lams1rus
+phirus(-1)*lams2rus*p2rus/mmrus(-1)
+phirus(-2)*lams3rus*p3rus/(mmrus(-1)*mmrus(-2))
+phirus(-3)*lams4rus*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+phirus(-4)*lams5rus*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+phirus(-5)*lams6rus*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+phirus(-6)*lams7rus*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+phirus(-7)*lams8rus*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7))
+(1-eduurus)*(1-phirus)*lamu1rus
+(1-phirus(-1))*lamu2rus*p2rus/mmrus(-1)
+(1-phirus(-2))*lamu3rus*p3rus/(mmrus(-1)*mmrus(-2))
+(1-phirus(-3))*lamu4rus*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+(1-phirus(-4))*lamu5rus*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+(1-phirus(-5))*lamu6rus*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+(1-phirus(-6))*lamu7rus*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+(1-phirus(-7))*lamu8rus*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)));

lablabsshchi=labchi/((1-eduschi)*phichi*lams1chi
+phichi(-1)*lams2chi*p2chi/mmchi(-1)
+phichi(-2)*lams3chi*p3chi/(mmchi(-1)*mmchi(-2))
+phichi(-3)*lams4chi*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+phichi(-4)*lams5chi*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+phichi(-5)*lams6chi*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+phichi(-6)*lams7chi*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+phichi(-7)*lams8chi*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7))
+(1-eduuchi)*(1-phichi)*lamu1chi
+(1-phichi(-1))*lamu2chi*p2chi/mmchi(-1)
+(1-phichi(-2))*lamu3chi*p3chi/(mmchi(-1)*mmchi(-2))
+(1-phichi(-3))*lamu4chi*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+(1-phichi(-4))*lamu5chi*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+(1-phichi(-5))*lamu6chi*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+(1-phichi(-6))*lamu7chi*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+(1-phichi(-7))*lamu8chi*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)));

lablabsshind=labind/((1-edusind)*phiind*lams1ind
+phiind(-1)*lams2ind*p2ind/mmind(-1)
+phiind(-2)*lams3ind*p3ind/(mmind(-1)*mmind(-2))
+phiind(-3)*lams4ind*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+phiind(-4)*lams5ind*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+phiind(-5)*lams6ind*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+phiind(-6)*lams7ind*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+phiind(-7)*lams8ind*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7))
+(1-eduuind)*(1-phiind)*lamu1ind
+(1-phiind(-1))*lamu2ind*p2ind/mmind(-1)
+(1-phiind(-2))*lamu3ind*p3ind/(mmind(-1)*mmind(-2))
+(1-phiind(-3))*lamu4ind*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+(1-phiind(-4))*lamu5ind*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+(1-phiind(-5))*lamu6ind*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+(1-phiind(-6))*lamu7ind*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+(1-phiind(-7))*lamu8ind*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)));

// ---- tx de prélèvement --- (labortax + constax) over gdp  ---------- //

taxesgdpadv=(tauadv*(wagddsadv*labsadv+wagdduadv*labdduadv)
+tcadv*(phiadv*xs1adv
+xs2adv*phiadv(-1)*p2adv/mmadv(-1)
+xs3adv*phiadv(-2)*p3adv/(mmadv(-1)*mmadv(-2))
+xs4adv*phiadv(-3)*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+xs5adv*phiadv(-4)*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+xs6adv*phiadv(-5)*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+xs7adv*phiadv(-6)*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+xs8adv*phiadv(-7)*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7)))
+tcadv*((1-phiadv)*xu1adv
+xu2adv*(1-phiadv(-1))*p2adv/mmadv(-1)
+xu3adv*(1-phiadv(-2))*p3adv/(mmadv(-1)*mmadv(-2))
+xu4adv*(1-phiadv(-3))*p4adv/(mmadv(-1)*mmadv(-2)*mmadv(-3))
+xu5adv*(1-phiadv(-4))*p5adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4))
+xu6adv*(1-phiadv(-5))*p6adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5))
+xu7adv*(1-phiadv(-6))*p7adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6))
+xu8adv*(1-phiadv(-7))*p8adv/(mmadv(-1)*mmadv(-2)*mmadv(-3)*mmadv(-4)*mmadv(-5)*mmadv(-6)*mmadv(-7))))/gdpadv;

taxesgdpnam=(taunam*(wagddsnam*labsnam+wagddunam*labddunam)
+tcnam*(phinam*xs1nam
+xs2nam*phinam(-1)*p2nam/mmnam(-1)
+xs3nam*phinam(-2)*p3nam/(mmnam(-1)*mmnam(-2))
+xs4nam*phinam(-3)*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+xs5nam*phinam(-4)*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+xs6nam*phinam(-5)*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+xs7nam*phinam(-6)*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+xs8nam*phinam(-7)*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7)))
+tcnam*((1-phinam)*xu1nam
+xu2nam*(1-phinam(-1))*p2nam/mmnam(-1)
+xu3nam*(1-phinam(-2))*p3nam/(mmnam(-1)*mmnam(-2))
+xu4nam*(1-phinam(-3))*p4nam/(mmnam(-1)*mmnam(-2)*mmnam(-3))
+xu5nam*(1-phinam(-4))*p5nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4))
+xu6nam*(1-phinam(-5))*p6nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5))
+xu7nam*(1-phinam(-6))*p7nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6))
+xu8nam*(1-phinam(-7))*p8nam/(mmnam(-1)*mmnam(-2)*mmnam(-3)*mmnam(-4)*mmnam(-5)*mmnam(-6)*mmnam(-7))))/gdpnam;

taxesgdpchi=(tauchi*(wagddschi*labschi+wagdduchi*labdduchi)
+tcchi*(phichi*xs1chi
+xs2chi*phichi(-1)*p2chi/mmchi(-1)
+xs3chi*phichi(-2)*p3chi/(mmchi(-1)*mmchi(-2))
+xs4chi*phichi(-3)*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+xs5chi*phichi(-4)*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+xs6chi*phichi(-5)*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+xs7chi*phichi(-6)*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+xs8chi*phichi(-7)*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7)))
+tcchi*((1-phichi)*xu1chi
+xu2chi*(1-phichi(-1))*p2chi/mmchi(-1)
+xu3chi*(1-phichi(-2))*p3chi/(mmchi(-1)*mmchi(-2))
+xu4chi*(1-phichi(-3))*p4chi/(mmchi(-1)*mmchi(-2)*mmchi(-3))
+xu5chi*(1-phichi(-4))*p5chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4))
+xu6chi*(1-phichi(-5))*p6chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5))
+xu7chi*(1-phichi(-6))*p7chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6))
+xu8chi*(1-phichi(-7))*p8chi/(mmchi(-1)*mmchi(-2)*mmchi(-3)*mmchi(-4)*mmchi(-5)*mmchi(-6)*mmchi(-7))))/gdpchi;

taxesgdpeas=(taueas*(wagddseas*labseas+wagddueas*labddueas)
+tceas*(phieas*xs1eas
+xs2eas*phieas(-1)*p2eas/mmeas(-1)
+xs3eas*phieas(-2)*p3eas/(mmeas(-1)*mmeas(-2))
+xs4eas*phieas(-3)*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+xs5eas*phieas(-4)*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+xs6eas*phieas(-5)*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+xs7eas*phieas(-6)*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+xs8eas*phieas(-7)*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7)))
+tceas*((1-phieas)*xu1eas
+xu2eas*(1-phieas(-1))*p2eas/mmeas(-1)
+xu3eas*(1-phieas(-2))*p3eas/(mmeas(-1)*mmeas(-2))
+xu4eas*(1-phieas(-3))*p4eas/(mmeas(-1)*mmeas(-2)*mmeas(-3))
+xu5eas*(1-phieas(-4))*p5eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4))
+xu6eas*(1-phieas(-5))*p6eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5))
+xu7eas*(1-phieas(-6))*p7eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6))
+xu8eas*(1-phieas(-7))*p8eas/(mmeas(-1)*mmeas(-2)*mmeas(-3)*mmeas(-4)*mmeas(-5)*mmeas(-6)*mmeas(-7))))/gdpeas;

taxesgdpind=(tauind*(wagddsind*labsind+wagdduind*labdduind)
+tcind*(phiind*xs1ind
+xs2ind*phiind(-1)*p2ind/mmind(-1)
+xs3ind*phiind(-2)*p3ind/(mmind(-1)*mmind(-2))
+xs4ind*phiind(-3)*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+xs5ind*phiind(-4)*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+xs6ind*phiind(-5)*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+xs7ind*phiind(-6)*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+xs8ind*phiind(-7)*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7)))
+tcind*((1-phiind)*xu1ind
+xu2ind*(1-phiind(-1))*p2ind/mmind(-1)
+xu3ind*(1-phiind(-2))*p3ind/(mmind(-1)*mmind(-2))
+xu4ind*(1-phiind(-3))*p4ind/(mmind(-1)*mmind(-2)*mmind(-3))
+xu5ind*(1-phiind(-4))*p5ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4))
+xu6ind*(1-phiind(-5))*p6ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5))
+xu7ind*(1-phiind(-6))*p7ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6))
+xu8ind*(1-phiind(-7))*p8ind/(mmind(-1)*mmind(-2)*mmind(-3)*mmind(-4)*mmind(-5)*mmind(-6)*mmind(-7))))/gdpind;

taxesgdpjap=(taujap*(wagddsjap*labsjap+wagddujap*labddujap)
+tcjap*(phijap*xs1jap
+xs2jap*phijap(-1)*p2jap/mmjap(-1)
+xs3jap*phijap(-2)*p3jap/(mmjap(-1)*mmjap(-2))
+xs4jap*phijap(-3)*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+xs5jap*phijap(-4)*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+xs6jap*phijap(-5)*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+xs7jap*phijap(-6)*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+xs8jap*phijap(-7)*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7)))
+tcjap*((1-phijap)*xu1jap
+xu2jap*(1-phijap(-1))*p2jap/mmjap(-1)
+xu3jap*(1-phijap(-2))*p3jap/(mmjap(-1)*mmjap(-2))
+xu4jap*(1-phijap(-3))*p4jap/(mmjap(-1)*mmjap(-2)*mmjap(-3))
+xu5jap*(1-phijap(-4))*p5jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4))
+xu6jap*(1-phijap(-5))*p6jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5))
+xu7jap*(1-phijap(-6))*p7jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6))
+xu8jap*(1-phijap(-7))*p8jap/(mmjap(-1)*mmjap(-2)*mmjap(-3)*mmjap(-4)*mmjap(-5)*mmjap(-6)*mmjap(-7))))/gdpjap;

taxesgdplac=(taulac*(wagddslac*labslac+wagddulac*labddulac)
+tclac*(philac*xs1lac
+xs2lac*philac(-1)*p2lac/mmlac(-1)
+xs3lac*philac(-2)*p3lac/(mmlac(-1)*mmlac(-2))
+xs4lac*philac(-3)*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+xs5lac*philac(-4)*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+xs6lac*philac(-5)*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+xs7lac*philac(-6)*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+xs8lac*philac(-7)*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7)))
+tclac*((1-philac)*xu1lac
+xu2lac*(1-philac(-1))*p2lac/mmlac(-1)
+xu3lac*(1-philac(-2))*p3lac/(mmlac(-1)*mmlac(-2))
+xu4lac*(1-philac(-3))*p4lac/(mmlac(-1)*mmlac(-2)*mmlac(-3))
+xu5lac*(1-philac(-4))*p5lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4))
+xu6lac*(1-philac(-5))*p6lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5))
+xu7lac*(1-philac(-6))*p7lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6))
+xu8lac*(1-philac(-7))*p8lac/(mmlac(-1)*mmlac(-2)*mmlac(-3)*mmlac(-4)*mmlac(-5)*mmlac(-6)*mmlac(-7))))/gdplac;

taxesgdpmen=(taumen*(wagddsmen*labsmen+wagddumen*labddumen)
+tcmen*(phimen*xs1men
+xs2men*phimen(-1)*p2men/mmmen(-1)
+xs3men*phimen(-2)*p3men/(mmmen(-1)*mmmen(-2))
+xs4men*phimen(-3)*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+xs5men*phimen(-4)*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+xs6men*phimen(-5)*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+xs7men*phimen(-6)*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+xs8men*phimen(-7)*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7)))
+tcmen*((1-phimen)*xu1men
+xu2men*(1-phimen(-1))*p2men/mmmen(-1)
+xu3men*(1-phimen(-2))*p3men/(mmmen(-1)*mmmen(-2))
+xu4men*(1-phimen(-3))*p4men/(mmmen(-1)*mmmen(-2)*mmmen(-3))
+xu5men*(1-phimen(-4))*p5men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4))
+xu6men*(1-phimen(-5))*p6men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5))
+xu7men*(1-phimen(-6))*p7men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6))
+xu8men*(1-phimen(-7))*p8men/(mmmen(-1)*mmmen(-2)*mmmen(-3)*mmmen(-4)*mmmen(-5)*mmmen(-6)*mmmen(-7))))/gdpmen;

taxesgdprus=(taurus*(wagddsrus*labsrus+wagddurus*labddurus)
+tcrus*(phirus*xs1rus
+xs2rus*phirus(-1)*p2rus/mmrus(-1)
+xs3rus*phirus(-2)*p3rus/(mmrus(-1)*mmrus(-2))
+xs4rus*phirus(-3)*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+xs5rus*phirus(-4)*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+xs6rus*phirus(-5)*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+xs7rus*phirus(-6)*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+xs8rus*phirus(-7)*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7)))
+tcrus*((1-phirus)*xu1rus
+xu2rus*(1-phirus(-1))*p2rus/mmrus(-1)
+xu3rus*(1-phirus(-2))*p3rus/(mmrus(-1)*mmrus(-2))
+xu4rus*(1-phirus(-3))*p4rus/(mmrus(-1)*mmrus(-2)*mmrus(-3))
+xu5rus*(1-phirus(-4))*p5rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4))
+xu6rus*(1-phirus(-5))*p6rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5))
+xu7rus*(1-phirus(-6))*p7rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6))
+xu8rus*(1-phirus(-7))*p8rus/(mmrus(-1)*mmrus(-2)*mmrus(-3)*mmrus(-4)*mmrus(-5)*mmrus(-6)*mmrus(-7))))/gdprus;

taxesgdpssa=(taussa*(wagddsssa*labsssa+wagddussa*labddussa)
+tcssa*(phissa*xs1ssa
+xs2ssa*phissa(-1)*p2ssa/mmssa(-1)
+xs3ssa*phissa(-2)*p3ssa/(mmssa(-1)*mmssa(-2))
+xs4ssa*phissa(-3)*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+xs5ssa*phissa(-4)*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+xs6ssa*phissa(-5)*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+xs7ssa*phissa(-6)*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+xs8ssa*phissa(-7)*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7)))
+tcssa*((1-phissa)*xu1ssa
+xu2ssa*(1-phissa(-1))*p2ssa/mmssa(-1)
+xu3ssa*(1-phissa(-2))*p3ssa/(mmssa(-1)*mmssa(-2))
+xu4ssa*(1-phissa(-3))*p4ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3))
+xu5ssa*(1-phissa(-4))*p5ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4))
+xu6ssa*(1-phissa(-5))*p6ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5))
+xu7ssa*(1-phissa(-6))*p7ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6))
+xu8ssa*(1-phissa(-7))*p8ssa/(mmssa(-1)*mmssa(-2)*mmssa(-3)*mmssa(-4)*mmssa(-5)*mmssa(-6)*mmssa(-7))))/gdpssa;

// -------------------------- growth gdp -------------------------------//

growthadv=(gdpadv/gdpadv(-1)-1);
growthnam=(gdpnam/gdpnam(-1)-1);
growthjap=(gdpjap/gdpjap(-1)-1);
growthssa=(gdpssa/gdpssa(-1)-1);
growthlac=(gdplac/gdplac(-1)-1);
growthmen=(gdpmen/gdpmen(-1)-1);
growthrus=(gdprus/gdprus(-1)-1);
growtheas=(gdpeas/gdpeas(-1)-1);
growthchi=(gdpchi/gdpchi(-1)-1);
growthind=(gdpind/gdpind(-1)-1);

// ------------------------- ownership ratio ---------------------------//

ownadv=weaadv/kadv;
ownnam=weanam/knam;
ownjap=weajap/kjap;
ownssa=weassa/kssa;
ownlac=wealac/klac;
ownmen=weamen/kmen;
ownrus=wearus/krus;
owneas=weaeas/keas;
ownchi=weachi/kchi;
ownind=weaind/kind;

ownsum=ownadv+ownnam+ownjap+ownssa+ownlac+ownmen+ownrus+owneas+ownchi+ownind;

// ----------------------- foreign assets fa ------------------------- //

faadv=weaadv-kadv;
fanam=weanam-knam;
fajap=weajap-kjap;
falac=wealac-klac;
faeas=weaeas-keas;
famen=weamen-kmen;
farus=wearus-krus;
fachi=weachi-kchi;
faind=weaind-kind;
fassa=weassa-kssa;

fasum=faadv+fanam+fajap+falac+faeas+famen+farus+fachi+faind+fassa;

// ----------------------- courrent account ca ----------------------- //

caadv=faadv-faadv(-1)/(ggnam*mmadv(-1));
canam=fanam-fanam(-1)/(ggnam*mmnam(-1));
cajap=fajap-fajap(-1)/(ggnam*mmjap(-1));
calac=falac-falac(-1)/(ggnam*mmlac(-1));
caeas=faeas-faeas(-1)/(ggnam*mmeas(-1));
camen=famen-famen(-1)/(ggnam*mmmen(-1));
carus=farus-farus(-1)/(ggnam*mmrus(-1));
cachi=fachi-fachi(-1)/(ggnam*mmchi(-1));
caind=faind-faind(-1)/(ggnam*mmind(-1));
cassa=fassa-fassa(-1)/(ggnam*mmssa(-1));

casum=caadv+canam+cajap+calac+caeas+camen+carus+cachi+caind+cassa;

// ------------------------- dette intérieure ---------------------------//

zzdebtadv=ddyadv*gdpadv;
zzdebtnam=ddynam*gdpnam;
zzdebtjap=ddyjap*gdpjap;
zzdebtlac=ddylac*gdplac;
zzdebteas=ddyeas*gdpeas;
zzdebtmen=ddymen*gdpmen;
zzdebtrus=ddyrus*gdprus;
zzdebtchi=ddychi*gdpchi;
zzdebtind=ddyind*gdpind;
zzdebtssa=ddyssa*gdpssa;

zzdebtsum=zzdebtadv+zzdebtnam+zzdebtjap+zzdebtlac+zzdebteas+zzdebtmen+zzdebtrus+zzdebtchi+zzdebtind+zzdebtssa;

// ------------------------- zzownership ratio ---------------------------//

zzown2adv=weaadv/(kadv+zzdebtadv);
zzown2nam=weanam/(knam+zzdebtnam);
zzown2jap=weajap/(kjap+zzdebtjap);
zzown2ssa=weassa/(kssa+zzdebtssa);
zzown2lac=wealac/(klac+zzdebtlac);
zzown2men=weamen/(kmen+zzdebtmen);
zzown2rus=wearus/(krus+zzdebtrus);
zzown2eas=weaeas/(keas+zzdebteas);
zzown2chi=weachi/(kchi+zzdebtchi);
zzown2ind=weaind/(kind+zzdebtind);

zzown2sum=zzown2adv+zzown2nam+zzown2jap+zzown2ssa+zzown2lac+zzown2men+zzown2rus+zzown2eas+zzown2chi+zzown2ind;

// ----------------------- foreign assets zzfa2 ------------------------- //

zzfa2adv=weaadv-kadv-zzdebtadv;
zzfa2nam=weanam-knam-zzdebtnam;
zzfa2jap=weajap-kjap-zzdebtjap;
zzfa2lac=wealac-klac-zzdebtlac;
zzfa2eas=weaeas-keas-zzdebteas;
zzfa2men=weamen-kmen-zzdebtmen;
zzfa2rus=wearus-krus-zzdebtrus;
zzfa2chi=weachi-kchi-zzdebtchi;
zzfa2ind=weaind-kind-zzdebtind;
zzfa2ssa=weassa-kssa-zzdebtssa;

zzfa2sum=zzfa2adv+zzfa2nam+zzfa2jap+zzfa2lac+zzfa2eas+zzfa2men+zzfa2rus+zzfa2chi+zzfa2ind+zzfa2ssa;

// ----------------------- courrent account zzca2 ----------------------- //

zzca2adv=zzfa2adv-zzfa2adv(-1)/(ggnam*mmadv(-1));
zzca2nam=zzfa2nam-zzfa2nam(-1)/(ggnam*mmnam(-1));
zzca2jap=zzfa2jap-zzfa2jap(-1)/(ggnam*mmjap(-1));
zzca2lac=zzfa2lac-zzfa2lac(-1)/(ggnam*mmlac(-1));
zzca2eas=zzfa2eas-zzfa2eas(-1)/(ggnam*mmeas(-1));
zzca2men=zzfa2men-zzfa2men(-1)/(ggnam*mmmen(-1));
zzca2rus=zzfa2rus-zzfa2rus(-1)/(ggnam*mmrus(-1));
zzca2chi=zzfa2chi-zzfa2chi(-1)/(ggnam*mmchi(-1));
zzca2ind=zzfa2ind-zzfa2ind(-1)/(ggnam*mmind(-1));
zzca2ssa=zzfa2ssa-zzfa2ssa(-1)/(ggnam*mmssa(-1));

zzca2sum=zzca2adv+zzca2nam+zzca2jap+zzca2lac+zzca2eas+zzca2men+zzca2rus+zzca2chi+zzca2ind+zzca2ssa;

// -------------------- new courrent account --------------------------- //

zzvardebtadv=zzdebtadv-zzdebtadv(-1);
zzvardebtnam=zzdebtnam-zzdebtnam(-1);
zzvardebtjap=zzdebtjap-zzdebtjap(-1);
zzvardebtlac=zzdebtlac-zzdebtlac(-1);
zzvardebteas=zzdebteas-zzdebteas(-1);
zzvardebtmen=zzdebtmen-zzdebtmen(-1);
zzvardebtrus=zzdebtrus-zzdebtrus(-1);
zzvardebtind=zzdebtind-zzdebtind(-1);
zzvardebtchi=zzdebtchi-zzdebtchi(-1);
zzvardebtssa=zzdebtssa-zzdebtssa(-1);

zzvarkadv=kadv-kadv(-1);
zzvarknam=knam-knam(-1);
zzvarkjap=kjap-kjap(-1);
zzvarklac=klac-klac(-1);
zzvarkeas=keas-keas(-1);
zzvarkmen=kmen-kmen(-1);
zzvarkrus=krus-krus(-1);
zzvarkind=kind-kind(-1);
zzvarkchi=kchi-kchi(-1);
zzvarkssa=kssa-kssa(-1);

zzvarweaadv=weaadv-weaadv(-1);
zzvarweanam=weanam-weanam(-1);
zzvarweajap=weajap-weajap(-1);
zzvarwealac=wealac-wealac(-1);
zzvarweaeas=weaeas-weaeas(-1);
zzvarweamen=weamen-weamen(-1);
zzvarwearus=wearus-wearus(-1);
zzvarweaind=weaind-weaind(-1);
zzvarweachi=weachi-weachi(-1);
zzvarweassa=weassa-weassa(-1);

zzznewcaadv=zzvarkadv+zzvardebtadv-zzvarweaadv;
zzznewcanam=zzvarknam+zzvardebtnam-zzvarweanam;
zzznewcajap=zzvarkjap+zzvardebtjap-zzvarweajap;
zzznewcalac=zzvarklac+zzvardebtlac-zzvarwealac;
zzznewcaeas=zzvarkeas+zzvardebteas-zzvarweaeas;
zzznewcamen=zzvarkmen+zzvardebtmen-zzvarweamen;
zzznewcarus=zzvarkrus+zzvardebtrus-zzvarwearus;
zzznewcaind=zzvarkind+zzvardebtind-zzvarweaind;
zzznewcachi=zzvarkchi+zzvardebtchi-zzvarweachi;
zzznewcassa=zzvarkssa+zzvardebtssa-zzvarweassa;

zzznewcasum=zzznewcaadv+zzznewcanam+zzznewcajap+zzznewcalac+zzznewcaeas+zzznewcamen+zzznewcarus+zzznewcachi+zzznewcaind+zzznewcassa;

// -------------------------- conspubgdp -------------------------------//

conspubgdpadv=conspubadv*gdpadv;
conspubgdpnam=conspubnam*gdpnam;
conspubgdpjap=conspubjap*gdpjap;
conspubgdpssa=conspubssa*gdpssa;
conspubgdplac=conspublac*gdplac;
conspubgdpmen=conspubmen*gdpmen;
conspubgdprus=conspubrus*gdprus;
conspubgdpeas=conspubeas*gdpeas;
conspubgdpchi=conspubchi*gdpchi;
conspubgdpind=conspubind*gdpind;

end;



initval;

bsadv           =       0.0879639;
bschi           =       0.0021782;
bseas           =       0.022133;
bsind           =       0.00620355;
bsjap           =       0.0605122;
bslac           =       0.0100509;
bsmen           =       0.0134604;
bsnam           =       0.0751243;
bsrus           =       0.0171303;
bsssa           =       0.000954535;

buadv           =       0.0494179;
buchi           =       0.0021782;
bueas           =       0.022133;
buind           =       0.00620355;
bujap           =       0.0296629;
bulac           =       0.0100509;
bumen           =       0.0134604;
bunam           =       0.0596225;
burus           =       0.0171303;
bussa           =       0.000954535;

gdpadv          =       0.993448;
gdpchi          =       0.0690066;
gdpeas          =       0.561693;
gdpind          =       0.094746;
gdpjap          =       0.737872;
gdplac          =       0.342757;
gdpmen          =       0.343658;
gdpnam          =       1.13712;
gdprus          =       0.376223;
gdpssa          =       0.125961;

govadv          =       0.013303;
govchi          =       0.0000203725;
goveas          =       0.000382988;
govind          =       0.0000269601;
govjap          =       0.0150679;
govlac          =       0.000151691;
govmen          =       0.000197734;
govnam          =       0.0144935;
govrus          =       0.000336704;
govssa          =    0.0000546259;

intrate         =       1.49773;

kadv            =       0.365185;
kchi            =       0.0253664;
keas            =       0.206475;
kind            =       0.034828;
kjap            =       0.271237;
klac            =       0.125995;
kmen            =       0.126326;
knam            =       0.417997;
krus            =       0.138297;
kssa            =       0.0463025;

tauadv    =             0.220401;
tauchi     =            0.0724983;
taueas      =           0.144628;
tauind       =          0.130217;
taujap        =         0.146443;
taulac        =         0.107896;
taumen        =         0.230856;
taunam        =         0.237564;
taurus        =         0.178388;
taussa        =         0.0960451;

weaadv          =       0.573545;
weachi          =       0.00873825;
weaeas          =       0.256306;
weaind          =       0.0138542;
weajap          =       0.244631;
wealac          =       0.101088;
weamen          =       0.0876558;
weanam          =       0.506575;
wearus          =       0.134789;
weassa          =       0.0199077;

xs1adv          =       0.000399832;
xs1chi          =       0.00488419;
xs1eas          =       0.0763439;
xs1ind          =       0.00531336;
xs1jap          =       0.252527;
xs1lac          =       0.0597624;
xs1men          =       0.0387268;
xs1nam          =       0.320226;
xs1rus          =       0.0467355;
xs1ssa          =       0.0103293;

xs2adv          =       0.000500379;
xs2chi          =       0.00611243;
xs2eas          =       0.0955425;
xs2ind          =       0.00664953;
xs2jap          =       0.316031;
xs2lac          =       0.0747912;
xs2men          =       0.0484656;
xs2nam          =       0.400755;
xs2rus          =       0.0584883;
xs2ssa          =       0.0129268;

xs3adv          =       0.000626212;
xs3chi          =       0.00764956;
xs3eas          =       0.119569;
xs3ind          =       0.00832172;
xs3jap          =       0.395505;
xs3lac          =       0.0935993;
xs3men          =       0.0606534;
xs3nam          =       0.501535;
xs3rus          =       0.0731966;
xs3ssa          =       0.0161776;

xs4adv          =       0.000783688;
xs4chi          =       0.00957322;
xs4eas          =       0.149637;
xs4ind          =       0.0104144;
xs4jap          =       0.494964;
xs4lac          =       0.117137;
xs4men          =       0.0759062;
xs4nam          =       0.627658;
xs4rus          =       0.0916037;
xs4ssa          =       0.0202459;

xs5adv          =       0.000980766;
xs5chi          =       0.0119806;
xs5eas          =       0.187268;
xs5ind          =       0.0130334;
xs5jap          =       0.619435;
xs5lac          =       0.146594;
xs5men          =       0.0949947;
xs5nam          =       0.785498;
xs5rus          =       0.11464;
xs5ssa          =       0.0253372;

xs6adv          =       0.0012274;
xs6chi          =       0.0149935;
xs6eas          =       0.234361;
xs6ind          =       0.0163109;
xs6jap          =       0.775207;
xs6lac          =       0.183459;
xs6men          =       0.118884;
xs6nam          =       0.983031;
xs6rus          =       0.143469;
xs6ssa          =       0.0317089;

xs7adv          =       0.00153606;
xs7chi          =       0.018764;
xs7eas          =       0.293296;
xs7ind          =       0.0204127;
xs7jap          =       0.970152;
xs7lac          =       0.229594;
xs7men          =       0.14878;
xs7nam          =       1.23024;
xs7rus          =       0.179547;
xs7ssa          =       0.0396828;

xs8adv          =       0.00192235;
xs8chi          =       0.0234826;
xs8eas          =       0.367053;
xs8ind          =       0.025546;
xs8jap          =       1.21412;
xs8lac          =       0.287331;
xs8men          =       0.186194;
xs8nam          =       1.53961;
xs8rus          =       0.224699;
xs8ssa          =       0.0496621;

xu1adv          =       0.145054;
xu1chi          =       0.00286073;
xu1eas          =       0.0439333;
xu1ind          =       0.00301439;
xu1jap          =       0.145247;
xu1lac          =       0.0361301;
xu1men          =       0.0223818;
xu1nam          =       0.183847;
xu1rus          =       0.0275213;
xu1ssa          =       0.00595664;

xu2adv          =       0.181531;
xu2chi          =       0.00358013;
xu2eas          =       0.0549815;
xu2ind          =       0.00377243;
xu2jap          =       0.181773;
xu2lac          =       0.0452159;
xu2men          =       0.0280103;
xu2nam          =       0.23008;
xu2rus          =       0.0344422;
xu2ssa          =       0.00745458;

xu3adv          =       0.227182;
xu3chi          =       0.00448044;
xu3eas          =       0.0688079;
xu3ind          =       0.0047211;
xu3jap          =       0.227484;
xu3lac          =       0.0565865;
xu3men          =       0.0350541;
xu3nam          =       0.28794;
xu3rus          =       0.0431035;
xu3ssa          =       0.00932922;

xu4adv          =       0.284312;
xu4chi          =       0.00560715;
xu4eas          =       0.0861113;
xu4ind          =       0.00590833;
xu4jap          =       0.284691;
xu4lac          =       0.0708166;
xu4men          =       0.0438694;
xu4nam          =       0.360349;
xu4rus          =       0.0539429;
xu4ssa          =       0.0116753;

xu5adv          =       0.355809;
xu5chi          =       0.00701721;
xu5eas          =       0.107766;
xu5ind          =       0.00739413;
xu5jap          =       0.356283;
xu5lac          =       0.0886251;
xu5men          =       0.0549014;
xu5nam          =       0.450968;
xu5rus          =       0.0675082;
xu5ssa          =       0.0146113;

xu6adv          =       0.445286;
xu6chi          =       0.00878186;
xu6eas          =       0.134867;
xu6ind          =       0.00925356;
xu6jap          =       0.445879;
xu6lac          =       0.110912;
xu6men          =       0.0687077;
xu6nam          =       0.564375;
xu6rus          =       0.0844848;
xu6ssa          =       0.0182857;

xu7adv          =       0.557265;
xu7chi          =       0.0109903;
xu7eas          =       0.168782;
xu7ind          =       0.0115806;
xu7jap          =       0.558006;
xu7lac          =       0.138804;
xu7men          =       0.0859859;
xu7nam          =       0.706301;
xu7rus          =       0.105731;
xu7ssa          =       0.0228841;

xu8adv          =       0.697402;
xu8chi          =       0.013754;
xu8eas          =       0.211227;
xu8ind          =       0.0144928;
xu8jap          =       0.698331;
xu8lac          =       0.173709;
xu8men          =       0.107609;
xu8nam          =       0.883918;
xu8rus          =       0.132319;
xu8ssa          =       0.0286388;

zs2adv       =          -0.000769199;
zs2chi       =          -0.0000276826;
zs2eas       =          -0.000238203;
zs2ind       =          -0.00000847775;
zs2jap       =          -0.00200381;
zs2lac       =          -0.000211421;
zs2men       =          -0.0000608569;
zs2nam       =          -0.00531673;
zs2rus       =          -0.000270973;
zs2ssa       =          -0.0000153023;
zs3adv       =          0.00127683;
zs3chi       =          0.0000356709;
zs3eas       =          0.00100211;
zs3ind       =          0.0000219374;
zs3jap       =          0.00106761;
zs3lac       =          0.000252602;
zs3men       =          0.000121134;
zs3nam       =          0.00648828;
zs3rus       =          0.000848925;
zs3ssa       =          0.0000224589;
zs4adv       =          0.00265069;
zs4chi       =          0.0000759496;
zs4eas       =          0.00200234;
zs4ind       =          0.0000451285;
zs4jap       =          0.00259837;
zs4lac       =          0.000530455;
zs4men       =          0.000248052;
zs4nam       =          0.013711;
zs4rus       =          0.00169616;
zs4ssa       =          0.0000465965;
zs5adv       =          0.00343248;
zs5chi       =          0.0000940641;
zs5eas       =          0.00273632;
zs5ind       =          0.0000585869;
zs5jap       =          0.00314945;
zs5lac       =          0.000665846;
zs5men       =          0.000317401;
zs5nam       =          0.0177472;
zs5rus       =          0.00229048;
zs5ssa       =          0.0000576693;
zs6adv       =          0.0029937;
zs6chi       =          0.0000647072;
zs6eas       =          0.00223273;
zs6ind       =          0.0000431499;
zs6jap       =          0.00240093;
zs6lac       =          0.000476397;
zs6men       =          0.000230167;
zs6nam       =          0.0150142;
zs6rus       =          0.00185827;
zs6ssa       =          0.0000389781;
zs7adv       =          0.00129282;
zs7chi       =          0.0000214387;
zs7eas       =          0.000956671;
zs7ind       =          0.0000147898;
zs7jap       =          0.000835275;
zs7lac       =          0.0001775;
zs7men       =          0.0000788313;
zs7nam       =          0.00636224;
zs7rus       =          0.000833786;
zs7ssa       =          0.0000120502;
zs8adv       =          0.000243118;
zs8chi       =          0.00000256661;
zs8eas       =          0.000167295;
zs8ind       =          0.00000189823;
zs8jap       =          0.000127513;
zs8lac       =          0.000027157;
zs8men       =          0.00000861939;
zs8nam       =          0.00141281;
zs8rus       =          0.000174784;
zs8ssa       =          0.00000125136;

zu2adv          =       0.065095;
zu2chi          =       0.00113188;
zu2eas          =       0.0242198;
zu2ind          =       0.00148454;
zu2jap          =       0.0260368;
zu2lac          =       0.0133506;
zu2men          =       0.0100839;
zu2nam          =       0.0420942;
zu2rus          =       0.0136097;
zu2ssa          =       0.00247341;
zu3adv          =       0.104583;
zu3chi          =       0.00180724;
zu3eas          =       0.0422828;
zu3ind          =       0.00260809;
zu3jap          =       0.0495755;
zu3lac          =       0.0204486;
zu3men          =       0.016935;
zu3nam          =       0.082379;
zu3rus          =       0.0227462;
zu3ssa          =       0.00407417;
zu4adv          =       0.12704;
zu4chi          =       0.00205922;
zu4eas          =       0.0547875;
zu4ind          =       0.00324124;
zu4jap          =       0.0577143;
zu4lac          =       0.0230174;
zu4men          =       0.0203709;
zu4nam          =       0.103604;
zu4rus          =       0.0285324;
zu4ssa          =       0.00472059;
zu5adv          =       0.131196;
zu5chi          =       0.00196138;
zu5eas          =       0.0603009;
zu5ind          =       0.0033146;
zu5jap          =       0.0553282;
zu5lac          =       0.0222218;
zu5men          =       0.0206023;
zu5nam          =       0.107669;
zu5rus          =       0.0304534;
zu5ssa          =       0.00465565;
zu6adv          =       0.0962605;
zu6chi          =       0.00124973;
zu6eas          =       0.0459305;
zu6ind          =       0.00226309;
zu6jap          =       0.034952;
zu6lac          =       0.0146768;
zu6men          =       0.0138615;
zu6nam          =       0.0764772;
zu6rus          =       0.022965;
zu6ssa          =       0.00293021;
zu7adv          =       0.0416822;
zu7chi          =       0.000417616;
zu7eas          =       0.0199463;
zu7ind          =       0.000795451;
zu7jap          =       0.0122058;
zu7lac          =       0.00550799;
zu7men          =       0.00481262;
zu7nam          =       0.0329505;
zu7rus          =       0.0104473;
zu7ssa          =       0.000909698;
zu8adv          =       0.00780175;
zu8chi          =       0.0000503712;
zu8eas          =       0.00352983;
zu8ind          =       0.000104291;
zu8jap          =       0.00185819;
zu8lac          =       0.00084801;
zu8men          =       0.000532558;
zu8nam          =       0.00738873;
zu8rus          =       0.00221635;
zu8ssa          =       0.0000948063;

// -------------------------------- open ------------------------------- //

ggnam=1.2;

// -------------------------- advanced countries -------------------- //

p2adv=  0.737909705     ;
p3adv=  0.583190133     ;
p4adv=  0.503607555     ;
p5adv=  0.378770303     ;
p6adv=  0.246313345     ;
p7adv=  0.106602647     ;
p8adv=  0.018188066     ;

edusadv=0.6;
eduuadv=0.2;

lams1adv=1;
lams2adv=1;
lams3adv=1;
lams4adv=1;
lams5adv=0.7;
lams6adv=0;
lams7adv=0;
lams8adv=0;

lamu1adv=1;
lamu2adv=1;
lamu3adv=1;
lamu4adv=1;
lamu5adv=0.5;
lamu6adv=0;
lamu7adv=0;
lamu8adv=0;

// -------------------------------- nam ------------------------------- //

p2nam=  0.636432622     ;
p3nam=  0.461119686     ;
p4nam=  0.390384311     ;
p5nam=  0.283128152     ;
p6nam=  0.170478555     ;
p7nam=  0.07064341      ;
p8nam=  0.014933077     ;

edusnam=0.6;
eduunam=0.2;

lams1nam=1;
lams2nam=1;
lams3nam=1;
lams4nam=1;
lams5nam=0.7;
lams6nam=0;
lams7nam=0;
lams8nam=0;

lamu1nam=1;
lamu2nam=1;
lamu3nam=1;
lamu4nam=1;
lamu5nam=0.5;
lamu6nam=0;
lamu7nam=0;
lamu8nam=0;

// -------------------------------- ssa ------------------------------- //

p2ssa=  0.902353166     ;
p3ssa=  0.738544114     ;
p4ssa=  0.498549334     ;
p5ssa=  0.338863948     ;
p6ssa=  0.18145258      ;
p7ssa=  0.05461183      ;
p8ssa=  0.00494107      ;

edusssa=0.6;
eduussa=0;

lams1ssa=1;
lams2ssa=1;
lams3ssa=1;
lams4ssa=1;
lams5ssa=0.5;
lams6ssa=0;
lams7ssa=0;
lams8ssa=0;

lamu1ssa=1;
lamu2ssa=1;
lamu3ssa=1;
lamu4ssa=1;
lamu5ssa=0.5;
lamu6ssa=0;
lamu7ssa=0;
lamu8ssa=0;

// ------------------------------- lac ------------------------------- //

p2lac=  0.750142354     ;
p3lac=  0.556601439     ;
p4lac=  0.395251537     ;
p5lac=  0.262137139     ;
p6lac=  0.147394629     ;
p7lac=  0.054555101     ;
p8lac=  0.007563576     ;

eduslac=0.6;
eduulac=0;

lams1lac=1;
lams2lac=1;
lams3lac=1;
lams4lac=1;
lams5lac=0.5;
lams6lac=0;
lams7lac=0;
lams8lac=0;

lamu1lac=1;
lamu2lac=1;
lamu3lac=1;
lamu4lac=1;
lamu5lac=0.5;
lamu6lac=0;
lamu7lac=0;
lamu8lac=0;

// ------------------------------- japan ------------------------------ //

p2jap=  0.6281131       ;
p3jap=  0.415373856     ;
p4jap=  0.289803583     ;
p5jap=  0.180504043     ;
p6jap=  0.093402302     ;
p7jap=  0.031027247     ;
p8jap=  0.004218908     ;

edusjap=0.6;
eduujap=0.2;

lams1jap=1;
lams2jap=1;
lams3jap=1;
lams4jap=1;
lams5jap=0.7;
lams6jap=0;
lams7jap=0;
lams8jap=0;

lamu1jap=1;
lamu2jap=1;
lamu3jap=1;
lamu4jap=1;
lamu5jap=0.5;
lamu6jap=0;
lamu7jap=0;
lamu8jap=0;

// ---------------------------- rusland ------------------------------- //

p2rus=  0.779655487     ;
p3rus=  0.651536269     ;
p4rus=  0.587062227     ;
p5rus=  0.457712116     ;
p6rus=  0.304849621     ;
p7rus=  0.139703701     ;
p8rus=  0.027761084     ;

edusrus=0.6;
eduurus=0;

lams1rus=1;
lams2rus=1;
lams3rus=1;
lams4rus=1;
lams5rus=0.5;
lams6rus=0;
lams7rus=0;
lams8rus=0;

lamu1rus=1;
lamu2rus=1;
lamu3rus=1;
lamu4rus=1;
lamu5rus=0.5;
lamu6rus=0;
lamu7rus=0;
lamu8rus=0;

// ------------------------------ mena -------------------------------- //

p2men=  0.883787741     ;
p3men=  0.745430149     ;
p4men=  0.568707292     ;
p5men=  0.419926284     ;
p6men=  0.248344451     ;
p7men=  0.084274197     ;
p8men=  0.007979153     ;

edusmen=0.6;
eduumen=0;

lams1men=1;
lams2men=1;
lams3men=1;
lams4men=1;
lams5men=0.5;
lams6men=0;
lams7men=0;
lams8men=0;

lamu1men=1;
lamu2men=1;
lamu3men=1;
lamu4men=1;
lamu5men=0.5;
lamu6men=0;
lamu7men=0;
lamu8men=0;

// ---------------------------- eastern europe ------------------------- //

p2eas=  0.847752891     ;
p3eas=  0.747812449     ;
p4eas=  0.685932238     ;
p5eas=  0.562651715     ;
p6eas=  0.385404591     ;
p7eas=  0.17026944      ;
p8eas=  0.027344764     ;

eduseas=0.6;
eduueas=0;

lams1eas=1;
lams2eas=1;
lams3eas=1;
lams4eas=1;
lams5eas=0.5;
lams6eas=0;
lams7eas=0;
lams8eas=0;

lamu1eas=1;
lamu2eas=1;
lamu3eas=1;
lamu4eas=1;
lamu5eas=0.5;
lamu6eas=0;
lamu7eas=0;
lamu8eas=0;

// ------------------------------ china ------------------------------- //

p2chi=  0.836132961     ;
p3chi=  0.683547181     ;
p4chi=  0.491653627     ;
p5chi=  0.308911372     ;
p6chi=  0.16753083      ;
p7chi=  0.054191479     ;
p8chi=  0.005700256     ;

eduschi=0.6;
eduuchi=0;

lams1chi=1;
lams2chi=1;
lams3chi=1;
lams4chi=1;
lams5chi=0.5;
lams6chi=0;
lams7chi=0;
lams8chi=0;

lamu1chi=1;
lamu2chi=1;
lamu3chi=1;
lamu4chi=1;
lamu5chi=0.5;
lamu6chi=0;
lamu7chi=0;
lamu8chi=0;

// ------------------------------ india ------------------------------- //

p2ind=  0.959325838     ;
p3ind=  0.877249138     ;
p4ind=  0.714826409     ;
p5ind=  0.555043112     ;
p6ind=  0.335815557     ;
p7ind=  0.110766231     ;
p8ind=  0.012399205     ;

edusind=0.6;
eduuind=0;

lams1ind=1;
lams2ind=1;
lams3ind=1;
lams4ind=1;
lams5ind=0.5;
lams6ind=0;
lams7ind=0;
lams8ind=0;

lamu1ind=1;
lamu2ind=1;
lamu3ind=1;
lamu4ind=1;
lamu5ind=0.5;
lamu6ind=0;
lamu7ind=0;
lamu8ind=0;

//------------------- initval for conspublique  ------------------------//

conspubnam=     0.170280977     ;
conspubadv=     0.142865046     ;
conspubeas=     0.104169516     ;
conspubmen=     0.167804208     ;
conspublac=     0.101727665     ;
conspubjap=     0.108820692     ;
conspubssa=     0.100729755     ;
conspubrus=     0.126599505     ;
conspubchi=     0.077531613     ;
conspubind=     0.080276755     ;

//----------------------- initval for debtratio  -----------------------//

ddynam= 0.042809925     ;
ddyadv= 0.044976211     ;
ddyeas= 0.002290146     ;
ddymen= 0.001932557     ;
ddylac= 0.001486455     ;
ddyjap= 0.068588        ;
ddyssa= 0.001456596     ;
ddyrus= 0.003005931     ;
ddychi= 0.00099159      ;
ddyind= 0.000955737     ;

// --------------------- initval for consumption tax ------------------//

tcnam=  0.111111111     ;
tcadv=  0.111111111     ;
tceas=  0.061111111     ;
tcmen=  0.061111111     ;
tclac=  0.061111111     ;
tcjap=  0.111111111     ;
tcssa=  0.061111111     ;
tcrus=  0.061111111     ;
tcchi=  0.061111111     ;
tcind=  0.061111111     ;

//-------------------- initval for country risk rating ----------------//

rankeas=        3.4     ;
rankmen=        3.952380952     ;
ranklac=        5.185185185     ;
rankssa=        6.404761905     ;
rankrus=        6.166666667     ;
rankchi=        3.181818182     ;
rankind=        4.888888889     ;

//-------------------- initval for risk premium pi ---------------------//

pimax=0.5;

pissa=rankssa/7*pimax;
pilac=ranklac/7*pimax;
pimen=rankmen/7*pimax;
pieas=rankeas/7*pimax;
pirus=rankrus/7*pimax;
pichi=rankchi/7*pimax;
piind=rankind/7*pimax;

//-------------------- initval for replacementratio  -------------------//

repnam= 0.2075  ;
repadv= 0.2125  ;
repeas= 0.2125  ;
repmen= 0.2     ;
replac= 0.125   ;
repjap= 0.1375  ;
repssa= 0.0375  ;
reprus= 0.225   ;
repchi= 0.15    ;
repind= 0.375   ;

//------------------- initval for skilled proportion  -------------------//

phinam= 0.093717733     ;
phiadv= 0.016504496     ;
phieas= 0.015494457     ;
phimen= 0.004347393     ;
philac= 0.009983981     ;
phijap= 0.03315712      ;
phissa= 0.002270809     ;
phirus= 0.02242942      ;
phichi= 0.006617984     ;
phiind= 0.002329718     ;

// ------------------------- initval for nn  -----------------------//

nnadvnam=       2.737660895     ;
nneasnam=       0.776607208     ;
nnmennam=       0.643436259     ;
nnlacnam=       0.628185755     ;
nnjapnam=       0.316280356     ;
nnssanam=       1.078307209     ;
nnrusnam=       1.263917419     ;
nnchinam=       4.749303854     ;
nnindnam=       4.054790781     ;

//-------------------------- initval for mm   --------------------------//

mmnam=1;
mmadv=1;
mmeas=1;
mmmen=1;
mmlac=1;
mmjap=1;
mmssa=1;
mmrus=1;
mmchi=1;
mmind=1;

//----------------------- initval for ratiogdp  ------------------------//

ratiogdpadv=0.739851327;
ratiogdpeas=0.337751537;
ratiogdpmen=0.231113615;
ratiogdplac=0.287509372;
ratiogdpjap=0.743360884;
ratiogdpssa=0.090156639;
ratiogdprus=0.253666004;
ratiogdpchi=0.051781166;
ratiogdpind=0.055246393;

//-------------------------- initval for psi  --------------------------//

psinam= 0.25    ;
psiadv= 0.25    ;
psieas= 0.25    ;
psimen= 0.25    ;
psilac= 0.25    ;
psijap= 0.25    ;
psissa= 0.25    ;
psirus= 0.25    ;
psichi= 0.25    ;
psiind= 0.25    ;

//-------------------------- initval for tfp  --------------------------//

tfpadv=         0.5486741       ;
tfpchi=         0.038072393     ;
tfpeas=         0.29634405      ;
tfpind=         0.044222983     ;
tfpjap=         0.57740563      ;
tfplac= 0.23485451      ;
tfpmen=         0.17735133      ;
tfprus=         0.25550939      ;
tfpssa=         0.072830822     ;


// -------------------------- non-walras ----------------------------- //

// ------------ initval --------- labs and labu ---------------------- //

labsadv         =       0.0410937;
labschi         =       0.0169804;
labseas         =       0.0459073;
labsind         =       0.00752248;
labsjap         =       0.0616604;
labslac         =       0.0222949;
labsmen         =       0.012207;
labsnam         =       0.195507;
labsrus         =       0.0593731;
labsssa         =       0.00615135;

labuadv         =       2.76765;
labuchi         =       3.14484;
labueas         =       3.50762;
labuind         =       3.82;
labujap         =       2.14982;
labulac         =       2.80478;
labumen         =       3.39307;
labunam         =       2.20181;
laburus         =       3.17428;
labussa         =       3.30136;

//------------------ initval for wags and wagu  ----------------------//

wagsadv       =         0.216341;
wagschi       =         0.0184062;
wagseas       =         0.122834;
wagsind       =         0.02097;
wagsjap       =         0.257108;
wagslac       =         0.0990729;
wagsmen       =         0.0788896;
wagsnam      =          0.359164;
wagsrus     =           0.0975531;
wagsssa    =            0.0350297;

waguadv   =             0.198091;
waguchi      =          0.0121864;
wagueas       =         0.0881035;
waguind       =         0.0138838;
wagujap       =         0.184412;
wagulac       =         0.0676766;
wagumen       =         0.0565839;
wagunam       =         0.257612;
wagurus       =         0.064588;
wagussa       =         0.0215359;

//-------------------- endval for  aa  ------------------------//

aaadv=          0.051225074     ;
aachi=          0.034984234     ;
aaeas=          0.059252538     ;
aaind=          0.017332098     ;
aajap=          0.099350937     ;
aalac=  0.044266179     ;
aamen=          0.024427738     ;
aanam=          0.19824889      ;
aarus=          0.080930529     ;
aassa=          0.01793136      ;


//-------------------- endval for   eta  ------------------------//

etaadv=         0.028799459     ;
etachi=         0.003354951     ;
etaeas=         0.005201302     ;
etaind=         0.006033906     ;
etajap=         0.001886245     ;
etalac= 0.00439495      ;
etamen=         0.007953791     ;
etanam=         0.00734676      ;
etarus=         0.004312647     ;
etassa=         0.011768173     ;


// ---------- initval ------------ wagddu -------------------------------//

wagdduadv=etaadv*wagsadv+(1-etaadv)*waguadv;
wagddunam=etanam*wagsnam+(1-etanam)*wagunam;
wagddussa=etassa*wagsssa+(1-etassa)*wagussa;
wagddulac=etalac*wagslac+(1-etalac)*wagulac;
wagddujap=etajap*wagsjap+(1-etajap)*wagujap;
wagddueas=etaeas*wagseas+(1-etaeas)*wagueas;
wagddumen=etamen*wagsmen+(1-etamen)*wagumen;
wagddurus=etarus*wagsrus+(1-etarus)*wagurus;
wagdduchi=etachi*wagschi+(1-etachi)*waguchi;
wagdduind=etaind*wagsind+(1-etaind)*waguind;

labdduadv     =         2.74665;
labdduchi     =         3.12882;
labddueas     =         3.48653;
labdduind     =         3.78506;
labddujap     =         2.1452;
labddulac     =         2.78787;
labddumen     =         3.36162;
labddunam     =         2.18494;
labddurus     =         3.154;
labddussa     =         3.23047;

wagddsadv=aaadv*(1-alpha)*kadv^alpha*tfpadv^(1-alpha)*labsadv^(sigma-1)
*(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^((1-alpha-sigma)/sigma);
wagddsnam=aanam*(1-alpha)*knam^alpha*labsnam^(sigma-1)
*(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^((1-alpha-sigma)/sigma);
wagddsssa=aassa*(1-alpha)*kssa^alpha*tfpssa^(1-alpha)*labsssa^(sigma-1)
*(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^((1-alpha-sigma)/sigma);
wagddslac=aalac*(1-alpha)*klac^alpha*tfplac^(1-alpha)*labslac^(sigma-1)
*(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^((1-alpha-sigma)/sigma);
wagddsjap=aajap*(1-alpha)*kjap^alpha*tfpjap^(1-alpha)*labsjap^(sigma-1)
*(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^((1-alpha-sigma)/sigma);
wagddseas=aaeas*(1-alpha)*keas^alpha*tfpeas^(1-alpha)*labseas^(sigma-1)
*(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^((1-alpha-sigma)/sigma);
wagddsmen=aamen*(1-alpha)*kmen^alpha*tfpmen^(1-alpha)*labsmen^(sigma-1)
*(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^((1-alpha-sigma)/sigma);
wagddsrus=aarus*(1-alpha)*krus^alpha*tfprus^(1-alpha)*labsrus^(sigma-1)
*(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^((1-alpha-sigma)/sigma);
wagddschi=aachi*(1-alpha)*kchi^alpha*tfpchi^(1-alpha)*labschi^(sigma-1)
*(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^((1-alpha-sigma)/sigma);
wagddsind=aaind*(1-alpha)*kind^alpha*tfpind^(1-alpha)*labsind^(sigma-1)
*(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^((1-alpha-sigma)/sigma);

// ------------- initval ------- unemeployement ----------------------- //

uradv=(labuadv-labdduadv)/labuadv;
urnam=(labunam-labddunam)/labunam;
urssa=(labussa-labddussa)/labussa;
urlac=(labulac-labddulac)/labulac;
urjap=(labujap-labddujap)/labujap;
ureas=(labueas-labddueas)/labueas;
urmen=(labumen-labddumen)/labumen;
urrus=(laburus-labddurus)/laburus;
urchi=(labuchi-labdduchi)/labuchi;
urind=(labuind-labdduind)/labuind;

// -------------------------------- h ------------------------------- //

hadv=wagsadv/waguadv;
hnam=wagsnam/wagunam;
hssa=wagsssa/wagussa;
hlac=wagslac/wagulac;
hjap=wagsjap/wagujap;
heas=wagseas/wagueas;
hmen=wagsmen/wagumen;
hrus=wagsrus/wagurus;
hchi=wagschi/waguchi;
hind=wagsind/waguind;

// ------------------------- labddu & wagddu ------------------------- //

thetaadv=wagddsadv/wagdduadv;
thetanam=wagddsnam/wagddunam;
thetassa=wagddsssa/wagddussa;
thetalac=wagddslac/wagddulac;
thetajap=wagddsjap/wagddujap;
thetaeas=wagddseas/wagddueas;
thetamen=wagddsmen/wagddumen;
thetarus=wagddsrus/wagddurus;
thetachi=wagddschi/wagdduchi;
thetaind=wagddsind/wagdduind;

// -------------------------------- ces ------------------------------- //

labadv=(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^(1/sigma);
labnam=(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^(1/sigma);
labssa=(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^(1/sigma);
lablac=(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^(1/sigma);
labjap=(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^(1/sigma);
labeas=(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^(1/sigma);
labmen=(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^(1/sigma);
labrus=(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^(1/sigma);
labchi=(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^(1/sigma);
labind=(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^(1/sigma);

// ---------- end --------- non-walras -------------------------------//

// ---------- fin --------- non-walras -------------------------------//

//-------------------- new variables ------ bengdp ---------------------//

bengdpadv       =       0.0281291;
bengdpchi       =       0.012054;
bengdpeas       =       0.0340587;
bengdpind       =       0.048223;
bengdpjap       =       0.00900521;
bengdplac       =       0.0099871;
bengdpmen       =       0.0215644;
bengdpnam       =       0.0210058;
bengdprus       =       0.0319259;
bengdpssa       =       0.0031103;

// ------------------ init ------ sup -------------------------------- //

supadv          =       8.63229;
supchi          =       14.5995;
supeas          =       6.59353;
supind          =       8.94687;
supjap          =       19.54;
suplac          =       14.1477;
supmen          =       10.6221;
supnam          =       10.8221;
suprus          =       7.35943;
supssa          =       14.4325;

// ------------------ init ------ trgdp ------------------------------- //


trgdpadv      =         0.0334561;
trgdpchi      =         0.0020779;
trgdpeas      =         0.00246767;
trgdpind      =         0.00235646;
trgdpjap      =         0.0299275;
trgdplac      =         0.00204177;
trgdpmen      =         0.00222994;
trgdpnam      =         0.0199405;
trgdprus      =         0.00234201;
trgdpssa      =         0.00211111;

// --------- initval ---------- rapport lab-lab sans h --------------- //

lablabsshadv  =         0.953313;
lablabsshchi  =         0.982495;
lablabssheas  =         0.958614;
lablabsshind  =         0.993544;
lablabsshjap  =         0.913683;
lablabsshlac  =         0.974443;
lablabsshmen  =         0.988267;
lablabsshnam  =         0.776518;
lablabsshrus  =         0.941908;
lablabsshssa  =         0.993889;

// --------- initval ---------- tx de prélèvement --------------- //

taxesgdpadv=(tauadv*(wagddsadv*labsadv+wagdduadv*labdduadv)
+tcadv*(phiadv*xs1adv
+xs2adv*phiadv*p2adv/mmadv
+xs3adv*phiadv*p3adv/(mmadv^2)
+xs4adv*phiadv*p4adv/(mmadv^3)
+xs5adv*phiadv*p5adv/(mmadv^4)
+xs6adv*phiadv*p6adv/(mmadv^5)
+xs7adv*phiadv*p7adv/(mmadv^6)
+xs8adv*phiadv*p8adv/(mmadv^7))
+tcadv*((1-phiadv)*xu1adv
+xu2adv*(1-phiadv)*p2adv/mmadv
+xu3adv*(1-phiadv)*p3adv/(mmadv^2)
+xu4adv*(1-phiadv)*p4adv/(mmadv^3)
+xu5adv*(1-phiadv)*p5adv/(mmadv^4)
+xu6adv*(1-phiadv)*p6adv/(mmadv^5)
+xu7adv*(1-phiadv)*p7adv/(mmadv^6)
+xu8adv*(1-phiadv)*p8adv/(mmadv^7)))/gdpadv;

taxesgdplac=(taulac*(wagddslac*labslac+wagddulac*labddulac)
+tclac*(philac*xs1lac
+xs2lac*philac*p2lac/mmlac
+xs3lac*philac*p3lac/(mmlac^2)
+xs4lac*philac*p4lac/(mmlac^3)
+xs5lac*philac*p5lac/(mmlac^4)
+xs6lac*philac*p6lac/(mmlac^5)
+xs7lac*philac*p7lac/(mmlac^6)
+xs8lac*philac*p8lac/(mmlac^7))
+tclac*((1-philac)*xu1lac
+xu2lac*(1-philac)*p2lac/mmlac
+xu3lac*(1-philac)*p3lac/(mmlac^2)
+xu4lac*(1-philac)*p4lac/(mmlac^3)
+xu5lac*(1-philac)*p5lac/(mmlac^4)
+xu6lac*(1-philac)*p6lac/(mmlac^5)
+xu7lac*(1-philac)*p7lac/(mmlac^6)
+xu8lac*(1-philac)*p8lac/(mmlac^7)))/gdplac;

taxesgdpjap=(taujap*(wagddsjap*labsjap+wagddujap*labddujap)
+tcjap*(phijap*xs1jap
+xs2jap*phijap*p2jap/mmjap
+xs3jap*phijap*p3jap/(mmjap^2)
+xs4jap*phijap*p4jap/(mmjap^3)
+xs5jap*phijap*p5jap/(mmjap^4)
+xs6jap*phijap*p6jap/(mmjap^5)
+xs7jap*phijap*p7jap/(mmjap^6)
+xs8jap*phijap*p8jap/(mmjap^7))
+tcjap*((1-phijap)*xu1jap
+xu2jap*(1-phijap)*p2jap/mmjap
+xu3jap*(1-phijap)*p3jap/(mmjap^2)
+xu4jap*(1-phijap)*p4jap/(mmjap^3)
+xu5jap*(1-phijap)*p5jap/(mmjap^4)
+xu6jap*(1-phijap)*p6jap/(mmjap^5)
+xu7jap*(1-phijap)*p7jap/(mmjap^6)
+xu8jap*(1-phijap)*p8jap/(mmjap^7)))/gdpjap;

taxesgdpeas=(taueas*(wagddseas*labseas+wagddueas*labddueas)
+tceas*(phieas*xs1eas
+xs2eas*phieas*p2eas/mmeas
+xs3eas*phieas*p3eas/(mmeas^2)
+xs4eas*phieas*p4eas/(mmeas^3)
+xs5eas*phieas*p5eas/(mmeas^4)
+xs6eas*phieas*p6eas/(mmeas^5)
+xs7eas*phieas*p7eas/(mmeas^6)
+xs8eas*phieas*p8eas/(mmeas^7))
+tceas*((1-phieas)*xu1eas
+xu2eas*(1-phieas)*p2eas/mmeas
+xu3eas*(1-phieas)*p3eas/(mmeas^2)
+xu4eas*(1-phieas)*p4eas/(mmeas^3)
+xu5eas*(1-phieas)*p5eas/(mmeas^4)
+xu6eas*(1-phieas)*p6eas/(mmeas^5)
+xu7eas*(1-phieas)*p7eas/(mmeas^6)
+xu8eas*(1-phieas)*p8eas/(mmeas^7)))/gdpeas;

taxesgdpnam=(taunam*(wagddsnam*labsnam+wagddunam*labddunam)
+tcnam*(phinam*xs1nam
+xs2nam*phinam*p2nam/mmnam
+xs3nam*phinam*p3nam/(mmnam^2)
+xs4nam*phinam*p4nam/(mmnam^3)
+xs5nam*phinam*p5nam/(mmnam^4)
+xs6nam*phinam*p6nam/(mmnam^5)
+xs7nam*phinam*p7nam/(mmnam^6)
+xs8nam*phinam*p8nam/(mmnam^7))
+tcnam*((1-phinam)*xu1nam
+xu2nam*(1-phinam)*p2nam/mmnam
+xu3nam*(1-phinam)*p3nam/(mmnam^2)
+xu4nam*(1-phinam)*p4nam/(mmnam^3)
+xu5nam*(1-phinam)*p5nam/(mmnam^4)
+xu6nam*(1-phinam)*p6nam/(mmnam^5)
+xu7nam*(1-phinam)*p7nam/(mmnam^6)
+xu8nam*(1-phinam)*p8nam/(mmnam^7)))/gdpnam;

taxesgdpmen=(taumen*(wagddsmen*labsmen+wagddumen*labddumen)
+tcmen*(phimen*xs1men
+xs2men*phimen*p2men/mmmen
+xs3men*phimen*p3men/(mmmen^2)
+xs4men*phimen*p4men/(mmmen^3)
+xs5men*phimen*p5men/(mmmen^4)
+xs6men*phimen*p6men/(mmmen^5)
+xs7men*phimen*p7men/(mmmen^6)
+xs8men*phimen*p8men/(mmmen^7))
+tcmen*((1-phimen)*xu1men
+xu2men*(1-phimen)*p2men/mmmen
+xu3men*(1-phimen)*p3men/(mmmen^2)
+xu4men*(1-phimen)*p4men/(mmmen^3)
+xu5men*(1-phimen)*p5men/(mmmen^4)
+xu6men*(1-phimen)*p6men/(mmmen^5)
+xu7men*(1-phimen)*p7men/(mmmen^6)
+xu8men*(1-phimen)*p8men/(mmmen^7)))/gdpmen;

taxesgdpchi=(tauchi*(wagddschi*labschi+wagdduchi*labdduchi)
+tcchi*(phichi*xs1chi
+xs2chi*phichi*p2chi/mmchi
+xs3chi*phichi*p3chi/(mmchi^2)
+xs4chi*phichi*p4chi/(mmchi^3)
+xs5chi*phichi*p5chi/(mmchi^4)
+xs6chi*phichi*p6chi/(mmchi^5)
+xs7chi*phichi*p7chi/(mmchi^6)
+xs8chi*phichi*p8chi/(mmchi^7))
+tcchi*((1-phichi)*xu1chi
+xu2chi*(1-phichi)*p2chi/mmchi
+xu3chi*(1-phichi)*p3chi/(mmchi^2)
+xu4chi*(1-phichi)*p4chi/(mmchi^3)
+xu5chi*(1-phichi)*p5chi/(mmchi^4)
+xu6chi*(1-phichi)*p6chi/(mmchi^5)
+xu7chi*(1-phichi)*p7chi/(mmchi^6)
+xu8chi*(1-phichi)*p8chi/(mmchi^7)))/gdpchi;

taxesgdpind=(tauind*(wagddsind*labsind+wagdduind*labdduind)
+tcind*(phiind*xs1ind
+xs2ind*phiind*p2ind/mmind
+xs3ind*phiind*p3ind/(mmind^2)
+xs4ind*phiind*p4ind/(mmind^3)
+xs5ind*phiind*p5ind/(mmind^4)
+xs6ind*phiind*p6ind/(mmind^5)
+xs7ind*phiind*p7ind/(mmind^6)
+xs8ind*phiind*p8ind/(mmind^7))
+tcind*((1-phiind)*xu1ind
+xu2ind*(1-phiind)*p2ind/mmind
+xu3ind*(1-phiind)*p3ind/(mmind^2)
+xu4ind*(1-phiind)*p4ind/(mmind^3)
+xu5ind*(1-phiind)*p5ind/(mmind^4)
+xu6ind*(1-phiind)*p6ind/(mmind^5)
+xu7ind*(1-phiind)*p7ind/(mmind^6)
+xu8ind*(1-phiind)*p8ind/(mmind^7)))/gdpind;

taxesgdprus=(taurus*(wagddsrus*labsrus+wagddurus*labddurus)
+tcrus*(phirus*xs1rus
+xs2rus*phirus*p2rus/mmrus
+xs3rus*phirus*p3rus/(mmrus^2)
+xs4rus*phirus*p4rus/(mmrus^3)
+xs5rus*phirus*p5rus/(mmrus^4)
+xs6rus*phirus*p6rus/(mmrus^5)
+xs7rus*phirus*p7rus/(mmrus^6)
+xs8rus*phirus*p8rus/(mmrus^7))
+tcrus*((1-phirus)*xu1rus
+xu2rus*(1-phirus)*p2rus/mmrus
+xu3rus*(1-phirus)*p3rus/(mmrus^2)
+xu4rus*(1-phirus)*p4rus/(mmrus^3)
+xu5rus*(1-phirus)*p5rus/(mmrus^4)
+xu6rus*(1-phirus)*p6rus/(mmrus^5)
+xu7rus*(1-phirus)*p7rus/(mmrus^6)
+xu8rus*(1-phirus)*p8rus/(mmrus^7)))/gdprus;

taxesgdpssa=(taussa*(wagddsssa*labsssa+wagddussa*labddussa)
+tcssa*(phissa*xs1ssa
+xs2ssa*phissa*p2ssa/mmssa
+xs3ssa*phissa*p3ssa/(mmssa^2)
+xs4ssa*phissa*p4ssa/(mmssa^3)
+xs5ssa*phissa*p5ssa/(mmssa^4)
+xs6ssa*phissa*p6ssa/(mmssa^5)
+xs7ssa*phissa*p7ssa/(mmssa^6)
+xs8ssa*phissa*p8ssa/(mmssa^7))
+tcssa*((1-phissa)*xu1ssa
+xu2ssa*(1-phissa)*p2ssa/mmssa
+xu3ssa*(1-phissa)*p3ssa/(mmssa^2)
+xu4ssa*(1-phissa)*p4ssa/(mmssa^3)
+xu5ssa*(1-phissa)*p5ssa/(mmssa^4)
+xu6ssa*(1-phissa)*p6ssa/(mmssa^5)
+xu7ssa*(1-phissa)*p7ssa/(mmssa^6)
+xu8ssa*(1-phissa)*p8ssa/(mmssa^7)))/gdpssa;

// --------- init ----------- growth gdp ------------------------------//

growthadv=(gdpadv/gdpadv-1);
growthnam=(gdpnam/gdpnam-1);
growthjap=(gdpjap/gdpjap-1);
growthssa=(gdpssa/gdpssa-1);
growthlac=(gdplac/gdplac-1);
growthmen=(gdpmen/gdpmen-1);
growthrus=(gdprus/gdprus-1);
growtheas=(gdpeas/gdpeas-1);
growthchi=(gdpchi/gdpchi-1);
growthind=(gdpind/gdpind-1);

// -------- init ---------- ownership ratio ---------------------------//

ownadv=weaadv/kadv;
ownnam=weanam/knam;
ownjap=weajap/kjap;
ownssa=weassa/kssa;
ownlac=wealac/klac;
ownmen=weamen/kmen;
ownrus=wearus/krus;
owneas=weaeas/keas;
ownchi=weachi/kchi;
ownind=weaind/kind;

ownsum=ownadv+ownnam+ownjap+ownssa+ownlac+ownmen+ownrus+owneas+ownchi+ownind;

// -------- init ---------- foreign assets fa ------------------------- //

faadv=weaadv-kadv;
fanam=weanam-knam;
fajap=weajap-kjap;
falac=wealac-klac;
faeas=weaeas-keas;
famen=weamen-kmen;
farus=wearus-krus;
fachi=weachi-kchi;
faind=weaind-kind;
fassa=weassa-kssa;

fasum=faadv+fanam+fajap+falac+faeas+famen+farus+fachi+faind+fassa;

// --------- init --------- courrent account ca ----------------------- //

caadv=faadv-faadv/(ggnam*mmadv);
canam=fanam-fanam/(ggnam*mmnam);
cajap=fajap-fajap/(ggnam*mmjap);
calac=falac-falac/(ggnam*mmlac);
caeas=faeas-faeas/(ggnam*mmeas);
camen=famen-famen/(ggnam*mmmen);
carus=farus-farus/(ggnam*mmrus);
cachi=fachi-fachi/(ggnam*mmchi);
caind=faind-faind/(ggnam*mmind);
cassa=fassa-fassa/(ggnam*mmssa);

casum=caadv+canam+cajap+calac+caeas+camen+carus+cachi+caind+cassa;

// ------------------------- dette intérieure ---------------------------//

zzdebtadv=ddyadv*gdpadv;
zzdebtnam=ddynam*gdpnam;
zzdebtjap=ddyjap*gdpjap;
zzdebtlac=ddylac*gdplac;
zzdebteas=ddyeas*gdpeas;
zzdebtmen=ddymen*gdpmen;
zzdebtrus=ddyrus*gdprus;
zzdebtchi=ddychi*gdpchi;
zzdebtind=ddyind*gdpind;
zzdebtssa=ddyssa*gdpssa;

zzdebtsum=zzdebtadv+zzdebtnam+zzdebtjap+zzdebtlac+zzdebteas+zzdebtmen+zzdebtrus+zzdebtchi+zzdebtind+zzdebtssa;

// ------------------------- zzownership ratio ---------------------------//

zzown2adv=weaadv/(kadv+zzdebtadv);
zzown2nam=weanam/(knam+zzdebtnam);
zzown2jap=weajap/(kjap+zzdebtjap);
zzown2ssa=weassa/(kssa+zzdebtssa);
zzown2lac=wealac/(klac+zzdebtlac);
zzown2men=weamen/(kmen+zzdebtmen);
zzown2rus=wearus/(krus+zzdebtrus);
zzown2eas=weaeas/(keas+zzdebteas);
zzown2chi=weachi/(kchi+zzdebtchi);
zzown2ind=weaind/(kind+zzdebtind);

zzown2sum=zzown2adv+zzown2nam+zzown2jap+zzown2ssa+zzown2lac+zzown2men+zzown2rus+zzown2eas+zzown2chi+zzown2ind;

// ----------------------- foreign assets zzfa2 ------------------------- //

zzfa2adv=weaadv-kadv-zzdebtadv;
zzfa2nam=weanam-knam-zzdebtnam;
zzfa2jap=weajap-kjap-zzdebtjap;
zzfa2lac=wealac-klac-zzdebtlac;
zzfa2eas=weaeas-keas-zzdebteas;
zzfa2men=weamen-kmen-zzdebtmen;
zzfa2rus=wearus-krus-zzdebtrus;
zzfa2chi=weachi-kchi-zzdebtchi;
zzfa2ind=weaind-kind-zzdebtind;
zzfa2ssa=weassa-kssa-zzdebtssa;

zzfa2sum=zzfa2adv+zzfa2nam+zzfa2jap+zzfa2lac+zzfa2eas+zzfa2men+zzfa2rus+zzfa2chi+zzfa2ind+zzfa2ssa;

// ----------------------- courrent account zzca2 ----------------------- //

zzca2adv=zzfa2adv-zzfa2adv/(ggnam*mmadv);
zzca2nam=zzfa2nam-zzfa2nam/(ggnam*mmnam);
zzca2jap=zzfa2jap-zzfa2jap/(ggnam*mmjap);
zzca2lac=zzfa2lac-zzfa2lac/(ggnam*mmlac);
zzca2eas=zzfa2eas-zzfa2eas/(ggnam*mmeas);
zzca2men=zzfa2men-zzfa2men/(ggnam*mmmen);
zzca2rus=zzfa2rus-zzfa2rus/(ggnam*mmrus);
zzca2chi=zzfa2chi-zzfa2chi/(ggnam*mmchi);
zzca2ind=zzfa2ind-zzfa2ind/(ggnam*mmind);
zzca2ssa=zzfa2ssa-zzfa2ssa/(ggnam*mmssa);

zzca2sum=zzca2adv+zzca2nam+zzca2jap+zzca2lac+zzca2eas+zzca2men+zzca2rus+zzca2chi+zzca2ind+zzca2ssa;

// -------------------- new courrent account --------------------------- //

zzvardebtadv=zzdebtadv-zzdebtadv;
zzvardebtnam=zzdebtnam-zzdebtnam;
zzvardebtjap=zzdebtjap-zzdebtjap;
zzvardebtlac=zzdebtlac-zzdebtlac;
zzvardebteas=zzdebteas-zzdebteas;
zzvardebtmen=zzdebtmen-zzdebtmen;
zzvardebtrus=zzdebtrus-zzdebtrus;
zzvardebtind=zzdebtind-zzdebtind;
zzvardebtchi=zzdebtchi-zzdebtchi;
zzvardebtssa=zzdebtssa-zzdebtssa;

zzvarkadv=kadv-kadv;
zzvarknam=knam-knam;
zzvarkjap=kjap-kjap;
zzvarklac=klac-klac;
zzvarkeas=keas-keas;
zzvarkmen=kmen-kmen;
zzvarkrus=krus-krus;
zzvarkind=kind-kind;
zzvarkchi=kchi-kchi;
zzvarkssa=kssa-kssa;

zzvarweaadv=weaadv-weaadv;
zzvarweanam=weanam-weanam;
zzvarweajap=weajap-weajap;
zzvarwealac=wealac-wealac;
zzvarweaeas=weaeas-weaeas;
zzvarweamen=weamen-weamen;
zzvarwearus=wearus-wearus;
zzvarweaind=weaind-weaind;
zzvarweachi=weachi-weachi;
zzvarweassa=weassa-weassa;

zzznewcaadv=zzvarkadv+zzvardebtadv-zzvarweaadv;
zzznewcanam=zzvarknam+zzvardebtnam-zzvarweanam;
zzznewcajap=zzvarkjap+zzvardebtjap-zzvarweajap;
zzznewcalac=zzvarklac+zzvardebtlac-zzvarwealac;
zzznewcaeas=zzvarkeas+zzvardebteas-zzvarweaeas;
zzznewcamen=zzvarkmen+zzvardebtmen-zzvarweamen;
zzznewcarus=zzvarkrus+zzvardebtrus-zzvarwearus;
zzznewcaind=zzvarkind+zzvardebtind-zzvarweaind;
zzznewcachi=zzvarkchi+zzvardebtchi-zzvarweachi;
zzznewcassa=zzvarkssa+zzvardebtssa-zzvarweassa;

zzznewcasum=zzznewcaadv+zzznewcanam+zzznewcajap+zzznewcalac+zzznewcaeas+zzznewcamen+zzznewcarus+zzznewcachi+zzznewcaind+zzznewcassa;

// -------------------------- conspubgdp -------------------------------//

conspubgdpadv=conspubadv*gdpadv;
conspubgdpnam=conspubnam*gdpnam;
conspubgdpjap=conspubjap*gdpjap;
conspubgdpssa=conspubssa*gdpssa;
conspubgdplac=conspublac*gdplac;
conspubgdpmen=conspubmen*gdpmen;
conspubgdprus=conspubrus*gdprus;
conspubgdpeas=conspubeas*gdpeas;
conspubgdpchi=conspubchi*gdpchi;
conspubgdpind=conspubind*gdpind;

end;

steady(solve_algo=3);
//check;



// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//

// ------------------------------ endval ------------------------------//




endval;

// ------------endo var---------------//


bsadv       =           0.348342;
bschi       =           0.00816161;
bseas       =           0.0386837;
bsind       =           0.00877731;
bsjap       =           0.258279;
bslac       =           0.0132808;
bsmen       =           0.0145322;
bsnam       =           0.218346;
bsrus       =           0.0156579;
bsssa       =           0.00051421;
buadv       =           0.165483;
buchi       =           0.00816161;
bueas       =           0.0386837;
buind       =           0.00877731;
bujap       =           0.112295;
bulac       =           0.0132808;
bumen       =           0.0145322;
bunam       =           0.173291;
burus       =           0.0156579;
bussa       =           0.00051421;

gdpadv      =           3.34309;
gdpchi      =           0.192905;
gdpeas      =           0.70728;
gdpind      =           0.0827861;
gdpjap      =           3.65576;
gdplac      =           0.389729;
gdpmen      =           0.271159;
gdpnam      =           4.37187;
gdprus      =           0.271172;
gdpssa      =           0.043046;

govadv      =           0.057936;
govchi      =           0.0000926278;
goveas      =           0.000646183;
govind      =           0.0000548723;
govjap      =           0.128242;
govlac      =           0.000418856;
govmen      =           0.000305629;
govnam      =           0.0698228;
govrus      =           0.000213442;
govssa      =           0.0000362923;

intrate     =           1.46185;

kadv        =           1.28006;
kchi        =           0.0738627;
keas        =           0.270815;
kind        =           0.0316985;
kjap        =           1.39978;
klac        =           0.149226;
kmen        =           0.103826;
knam        =           1.67397;
krus        =           0.103831;
kssa        =           0.0164821;

tauadv      =           0.570843;
tauchi      =           0.268101;
taueas      =           0.366485;
tauind      =           0.403354;
taujap      =           0.466199;
taulac      =           0.284974;
taumen      =           0.34579;
taunam      =           0.36596;
taurus      =           0.336702;
taussa      =           0.204874;

weaadv      =           0.909112;
weachi      =           0.117435;
weaeas      =           0.320515;
weaind      =           0.0249183;
weajap      =           1.66902;
wealac      =           0.243415;
weamen      =           0.128726;
weanam      =           2.07816;
wearus      =           0.111155;
weassa      =           0.0230922;

xs1adv      =           0.188576;
xs1chi      =           0.0189034;
xs1eas      =           0.0565248;
xs1ind      =           0.00741483;
xs1jap      =           0.219451;
xs1lac      =           0.0356519;
xs1men      =           0.0232345;
xs1nam      =           0.245332;
xs1rus      =           0.0238094;
xs1ssa      =           0.00573089;
xs2adv      =           0.229726;
xs2chi      =           0.0230283;
xs2eas      =           0.0688591;
xs2ind      =           0.00903283;
xs2jap      =           0.267337;
xs2lac      =           0.0434315;
xs2men      =           0.0283045;
xs2nam      =           0.298866;
xs2rus      =           0.0290049;
xs2ssa      =           0.00698143;
xs3adv      =           0.279854;
xs3chi      =           0.0280533;
xs3eas      =           0.0838849;
xs3ind      =           0.0110039;
xs3jap      =           0.325673;
xs3lac      =           0.0529087;
xs3men      =           0.0344808;
xs3nam      =           0.364082;
xs3rus      =           0.035334;
xs3ssa      =           0.00850485;
xs4adv      =           0.340921;
xs4chi      =           0.0341749;
xs4eas      =           0.102189;
xs4ind      =           0.013405;
xs4jap      =           0.396738;
xs4lac      =           0.0644539;
xs4men      =           0.0420048;
xs4nam      =           0.443528;
xs4rus      =           0.0430442;
xs4ssa      =           0.0103607;
xs5adv      =           0.415314;
xs5chi      =           0.0416321;
xs5eas      =           0.124488;
xs5ind      =           0.0163301;
xs5jap      =           0.48331;
xs5lac      =           0.0785184;
xs5men      =           0.0511707;
xs5nam      =           0.54031;
xs5rus      =           0.0524369;
xs5ssa      =           0.0126215;
xs6adv      =           0.505939;
xs6chi      =           0.0507167;
xs6eas      =           0.151653;
xs6ind      =           0.0198935;
xs6jap      =           0.588774;
xs6lac      =           0.0956519;
xs6men      =           0.0623366;
xs6nam      =           0.658211;
xs6rus      =           0.0638792;
xs6ssa      =           0.0153756;
xs7adv      =           0.61634;
xs7chi      =           0.0617836;
xs7eas      =           0.184745;
xs7ind      =           0.0242345;
xs7jap      =           0.71725;
xs7lac      =           0.116524;
xs7men      =           0.0759391;
xs7nam      =           0.80184;
xs7rus      =           0.0778182;
xs7ssa      =           0.0187308;
xs8adv      =           0.750832;
xs8chi      =           0.0752654;
xs8eas      =           0.225058;
xs8ind      =           0.0295227;
xs8jap      =           0.873761;
xs8lac      =           0.141951;
xs8men      =           0.0925098;
xs8nam      =           0.976809;
xs8rus      =           0.094799;
xs8ssa      =           0.022818;

xu1adv      =           0.0977925;
xu1chi      =           0.0106559;
xu1eas      =           0.0325568;
xu1ind      =           0.00443707;
xu1jap      =           0.111348;
xu1lac      =           0.0200241;
xu1men      =           0.0133057;
xu1nam      =           0.136053;
xu1rus      =           0.0136576;
xu1ssa      =           0.00316713;
xu2adv      =           0.119132;
xu2chi      =           0.0129811;
xu2eas      =           0.039661;
xu2ind      =           0.00540528;
xu2jap      =           0.135645;
xu2lac      =           0.0243935;
xu2men      =           0.0162091;
xu2nam      =           0.165741;
xu2rus      =           0.0166378;
xu2ssa      =           0.00385823;
xu3adv      =           0.145128;
xu3chi      =           0.0158137;
xu3eas      =           0.0483155;
xu3ind      =           0.00658476;
xu3jap      =           0.165245;
xu3lac      =           0.0297164;
xu3men      =           0.0197461;
xu3nam      =           0.201907;
xu3rus      =           0.0202683;
xu3ssa      =           0.00470013;
xu4adv      =           0.176796;
xu4chi      =           0.0192644;
xu4eas      =           0.0588584;
xu4ind      =           0.00802162;
xu4jap      =           0.201303;
xu4lac      =           0.0362008;
xu4men      =           0.0240549;
xu4nam      =           0.245965;
xu4rus      =           0.0246911;
xu4ssa      =           0.00572575;
xu5adv      =           0.215374;
xu5chi      =           0.0234681;
xu5eas      =           0.0717018;
xu5ind      =           0.00977202;
xu5jap      =           0.245229;
xu5lac      =           0.0441002;
xu5men      =           0.0293039;
xu5nam      =           0.299637;
xu5rus      =           0.0300789;
xu5ssa      =           0.00697516;
xu6adv      =           0.262371;
xu6chi      =           0.028589;
xu6eas      =           0.0873479;
xu6ind      =           0.0119044;
xu6jap      =           0.29874;
xu6lac      =           0.0537233;
xu6men      =           0.0356983;
xu6nam      =           0.365021;
xu6rus      =           0.0366424;
xu6ssa      =           0.00849721;
xu7adv      =           0.319623;
xu7chi      =           0.0348274;
xu7eas      =           0.106408;
xu7ind      =           0.014502;
xu7jap      =           0.363928;
xu7lac      =           0.0654463;
xu7men      =           0.043488;
xu7nam      =           0.444672;
xu7rus      =           0.0446382;
xu7ssa      =           0.0103514;
xu8adv      =           0.389368;
xu8chi      =           0.0424271;
xu8eas      =           0.129627;
xu8ind      =           0.0176665;
xu8jap      =           0.443341;
xu8lac      =           0.0797273;
xu8men      =           0.0529775;
xu8nam      =           0.541704;
xu8rus      =           0.0543787;
xu8ssa      =           0.0126102;

zs2adv      =           -0.00932546;
zs2chi      =           -0.000118167;
zs2eas      =           -0.000993383;
zs2ind      =           -0.000109752;
zs2jap      =           -0.00683773;
zs2lac      =           -0.000284554;
zs2men      =           -0.000298793;
zs2nam      =           -0.00970037;
zs2rus      =           -0.000755272;
zs2ssa      =           -0.0000403761;
zs3adv      =           0.0270161;
zs3chi      =           0.0016419;
zs3eas      =           0.0087632;
zs3ind      =           0.000410655;
zs3jap      =           0.0585838;
zs3lac      =           0.00482788;
zs3men      =           0.00272908;
zs3nam      =           0.11799;
zs3rus      =           0.00375978;
zs3ssa      =           0.000166243;
zs4adv      =           0.0581496;
zs4chi      =           0.00331119;
zs4eas      =           0.0178318;
zs4ind      =           0.000863409;
zs4jap      =           0.119815;
zs4lac      =           0.0097086;
zs4men      =           0.00554865;
zs4nam      =           0.238217;
zs4rus      =           0.00779536;
zs4ssa      =           0.000347761;
zs5adv      =           0.0808529;
zs5chi      =           0.00482829;
zs5eas      =           0.0258694;
zs5ind      =           0.00121377;
zs5jap      =           0.174414;
zs5lac      =           0.0141308;
zs5men      =           0.00806328;
zs5nam      =           0.349621;
zs5rus      =           0.0109728;
zs5ssa      =           0.000478034;
zs6adv      =           0.0877184;
zs6chi      =           0.0046579;
zs6eas      =           0.0248172;
zs6ind      =           0.00114184;
zs6jap      =           0.194938;
zs6lac      =           0.0138119;
zs6men      =           0.00768292;
zs6nam      =           0.384453;
zs6rus      =           0.0100367;
zs6ssa      =           0.000406518;
zs7adv      =           0.0676125;
zs7chi      =           0.00280069;
zs7eas      =           0.0147763;
zs7ind      =           0.000677689;
zs7jap      =           0.145891;
zs7lac      =           0.00880478;
zs7men      =           0.00443579;
zs7nam      =           0.256298;
zs7rus      =           0.00570608;
zs7ssa      =           0.000196048;
zs8adv      =           0.0321169;
zs8chi      =           0.000973923;
zs8eas      =           0.00500638;
zs8ind      =           0.000219513;
zs8jap      =           0.0732824;
zs8lac      =           0.00344264;
zs8men      =           0.00136606;
zs8nam      =           0.103887;
zs8rus      =           0.00182655;
zs8ssa      =           0.0000477022;

zu2adv      =           0.0276745;
zu2chi      =           0.00697915;
zu2eas      =           0.0168485;
zu2ind      =           0.00192522;
zu2jap      =           0.0407056;
zu2lac      =           0.0127909;
zu2men      =           0.00745783;
zu2nam      =           0.0312747;
zu2rus      =           0.0059286;
zu2ssa      =           0.00181839;
zu3adv      =           0.0663281;
zu3chi      =           0.0132624;
zu3eas      =           0.0313889;
zu3ind      =           0.00336202;
zu3jap      =           0.097917;
zu3lac      =           0.0244314;
zu3men      =           0.013941;
zu3nam      =           0.0747455;
zu3rus      =           0.0108041;
zu3ssa      =           0.00334136;
zu4adv      =           0.0983058;
zu4chi      =           0.0186499;
zu4eas      =           0.0430135;
zu4ind      =           0.00419973;
zu4jap      =           0.151047;
zu4lac      =           0.0345763;
zu4men      =           0.0191905;
zu4nam      =           0.112646;
zu4rus      =           0.0143709;
zu4ssa      =           0.00449726;
zu5adv      =           0.118732;
zu5chi      =           0.0226366;
zu5eas      =           0.0501596;
zu5ind      =           0.00415828;
zu5jap      =           0.196773;
zu5lac      =           0.042352;
zu5men      =           0.0225353;
zu5nam      =           0.141869;
zu5rus      =           0.0160961;
zu5ssa      =           0.00517425;
zu6adv      =           0.119711;
zu6chi      =           0.0206435;
zu6eas      =           0.0451531;
zu6ind      =           0.00357049;
zu6jap      =           0.199297;
zu6lac      =           0.0392683;
zu6men      =           0.0201576;
zu6nam      =           0.13906;
zu6rus      =           0.0137628;
zu6ssa      =           0.00415519;
zu7adv      =           0.0917906;
zu7chi      =           0.0126792;
zu7eas      =           0.0280445;
zu7ind      =           0.00242234;
zu7jap      =           0.149251;
zu7lac      =           0.0254583;
zu7men      =           0.0120783;
zu7nam      =           0.0971875;
zu7rus      =           0.00815216;
zu7ssa      =           0.0020125;
zu8adv      =           0.0424288;
zu8chi      =           0.00448912;
zu8eas      =           0.00983556;
zu8ind      =           0.000863094;
zu8jap      =           0.0739413;
zu8lac      =           0.0100958;
zu8men      =           0.003838;
zu8nam      =           0.0406161;
zu8rus      =           0.00269801;
zu8ssa      =           0.000491313;

// -------------------------------- open ------------------------------- //

ggnam=1.2;

// -------------------------- advanced countries -------------------- //

p2adv=  0.986289852     ;
p3adv=  0.955148399     ;
p4adv=  1       ;
p5adv=  1       ;
p6adv=  0.976430451     ;
p7adv=  0.75710297      ;
p8adv=  0.379784022     ;

edusadv=0.6;
eduuadv=0.2;

lams1adv=1;
lams2adv=1;
lams3adv=1;
lams4adv=1;
lams5adv=0.8;
lams6adv=0;
lams7adv=0;
lams8adv=0;

lamu1adv=1;
lamu2adv=1;
lamu3adv=1;
lamu4adv=1;
lamu5adv=0.6;
lamu6adv=0;
lamu7adv=0;
lamu8adv=0;

// -------------------------------- nam ------------------------------- //

p2nam=  0.986289852     ;
p3nam=  0.955148399     ;
p4nam=  1       ;
p5nam=  1       ;
p6nam=  0.98579086      ;
p7nam=  0.742740395     ;
p8nam=  0.344004851     ;

edusnam=0.6;
eduunam=0.2;

lams1nam=1;
lams2nam=1;
lams3nam=1;
lams4nam=1;
lams5nam=0.8;
lams6nam=0;
lams7nam=0;
lams8nam=0;

lamu1nam=1;
lamu2nam=1;
lamu3nam=1;
lamu4nam=1;
lamu5nam=0.6;
lamu6nam=0;
lamu7nam=0;
lamu8nam=0;

// -------------------------------- ssa ------------------------------- //

p2ssa=  0.986289852     ;
p3ssa=  0.955148399     ;
p4ssa=  0.759653392     ;
p5ssa=  0.634364253     ;
p6ssa=  0.463109632     ;
p7ssa=  0.241577047     ;
p8ssa=  0.059963112     ;

edusssa=0.6;
eduussa=0;

lams1ssa=1;
lams2ssa=1;
lams3ssa=1;
lams4ssa=1;
lams5ssa=0.5;
lams6ssa=0;
lams7ssa=0;
lams8ssa=0;

lamu1ssa=1;
lamu2ssa=1;
lamu3ssa=1;
lamu4ssa=1;
lamu5ssa=0.5;
lamu6ssa=0;
lamu7ssa=0;
lamu8ssa=0;

// ------------------------------- lac ------------------------------- //

p2lac=  0.986289852     ;
p3lac=  0.955148399     ;
p4lac=  0.918592142     ;
p5lac=  0.860181373     ;
p6lac=  0.743582371     ;
p7lac=  0.539856069     ;
p8lac=  0.24944583      ;

eduslac=0.6;
eduulac=0;

lams1lac=1;
lams2lac=1;
lams3lac=1;
lams4lac=1;
lams5lac=0.5;
lams6lac=0;
lams7lac=0;
lams8lac=0;

lamu1lac=1;
lamu2lac=1;
lamu3lac=1;
lamu4lac=1;
lamu5lac=0.5;
lamu6lac=0;
lamu7lac=0;
lamu8lac=0;

// ------------------------------- japan ------------------------------ //

p2jap=  0.986289852     ;
p3jap=  0.955148399     ;
p4jap=  1       ;
p5jap=  0.99240326      ;
p6jap=  0.926086957     ;
p7jap=  0.773898233     ;
p8jap=  0.483304477     ;

edusjap=0.6;
eduujap=0.2;

lams1jap=1;
lams2jap=1;
lams3jap=1;
lams4jap=1;
lams5jap=0.8;
lams6jap=0;
lams7jap=0;
lams8jap=0;

lamu1jap=1;
lamu2jap=1;
lamu3jap=1;
lamu4jap=1;
lamu5jap=0.6;
lamu6jap=0;
lamu7jap=0;
lamu8jap=0;

// ---------------------------- rusland ------------------------------- //

p2rus=  0.986289852     ;
p3rus=  0.955148399     ;
p4rus=  0.870982898     ;
p5rus=  0.784798781     ;
p6rus=  0.615551778     ;
p7rus=  0.376706946     ;
p8rus=  0.128615306     ;

edusrus=0.6;
eduurus=0;

lams1rus=1;
lams2rus=1;
lams3rus=1;
lams4rus=1;
lams5rus=0.5;
lams6rus=0;
lams7rus=0;
lams8rus=0;

lamu1rus=1;
lamu2rus=1;
lamu3rus=1;
lamu4rus=1;
lamu5rus=0.5;
lamu6rus=0;
lamu7rus=0;
lamu8rus=0;

// ------------------------------ mena -------------------------------- //

p2men=  0.986289852     ;
p3men=  0.955148399     ;
p4men=  0.965180767     ;
p5men=  0.909809812     ;
p6men=  0.782296592     ;
p7men=  0.501229312     ;
p8men=  0.162832053     ;

edusmen=0.6;
eduumen=0;

lams1men=1;
lams2men=1;
lams3men=1;
lams4men=1;
lams5men=0.5;
lams6men=0;
lams7men=0;
lams8men=0;

lamu1men=1;
lamu2men=1;
lamu3men=1;
lamu4men=1;
lamu5men=0.5;
lamu6men=0;
lamu7men=0;
lamu8men=0;

// ---------------------------- eastern europe ------------------------- //

p2eas=  0.986289852     ;
p3eas=  0.955148399     ;
p4eas=  0.959572785     ;
p5eas=  0.90729252      ;
p6eas=  0.778366337     ;
p7eas=  0.508173052     ;
p8eas=  0.187129118     ;

eduseas=0.6;
eduueas=0;

lams1eas=1;
lams2eas=1;
lams3eas=1;
lams4eas=1;
lams5eas=0.5;
lams6eas=0;
lams7eas=0;
lams8eas=0;

lamu1eas=1;
lamu2eas=1;
lamu3eas=1;
lamu4eas=1;
lamu5eas=0.5;
lamu6eas=0;
lamu7eas=0;
lamu8eas=0;

// ------------------------------ china ------------------------------- //

p2chi=  0.986289852     ;
p3chi=  0.955148399     ;
p4chi=  0.943012567     ;
p5chi=  0.891094232     ;
p6chi=  0.772421699     ;
p7chi=  0.521311169     ;
p8chi=  0.202793499     ;

eduschi=0.6;
eduuchi=0;

lams1chi=1;
lams2chi=1;
lams3chi=1;
lams4chi=1;
lams5chi=0.5;
lams6chi=0;
lams7chi=0;
lams8chi=0;

lamu1chi=1;
lamu2chi=1;
lamu3chi=1;
lamu4chi=1;
lamu5chi=0.5;
lamu6chi=0;
lamu7chi=0;
lamu8chi=0;

// ------------------------------ india ------------------------------- //

p2ind=  0.986289852     ;
p3ind=  0.955148399     ;
p4ind=  0.941906608     ;
p5ind=  0.873350269     ;
p6ind=  0.717620394     ;
p7ind=  0.442768491     ;
p8ind=  0.146184372     ;

edusind=0.6;
eduuind=0;

lams1ind=1;
lams2ind=1;
lams3ind=1;
lams4ind=1;
lams5ind=0.5;
lams6ind=0;
lams7ind=0;
lams8ind=0;

lamu1ind=1;
lamu2ind=1;
lamu3ind=1;
lamu4ind=1;
lamu5ind=0.5;
lamu6ind=0;
lamu7ind=0;
lamu8ind=0;

//------------------- endval for conspublique  ------------------------//

conspubnam=     0.149242875     ;
conspubadv=     0.171279136     ;
conspubeas=     0.164212094     ;
conspubmen=     0.154196071     ;
conspublac=     0.150334205     ;
conspubjap=     0.164831366     ;
conspubssa=     0.150356571     ;
conspubrus=     0.162501299     ;
conspubchi=     0.125176012     ;
conspubind=     0.110200901     ;

//------------------------- endval for debtratio  -----------------------//

ddynam= 0.060992156     ;
ddyadv= 0.066182791     ;
ddyeas= 0.003489057     ;
ddymen= 0.004304415     ;
ddylac= 0.004104364     ;
ddyjap= 0.13396642      ;
ddyssa= 0.003219776     ;
ddyrus= 0.003005931     ;
ddychi= 0.001833757     ;
ddyind= 0.002531276     ;

//------------------- endval for replacement ratio  --------------------//

repnam= 0.36    ;
repadv= 0.37    ;
repeas= 0.425   ;
repmen= 0.4     ;
replac= 0.25    ;
repjap= 0.22    ;
repssa= 0.075   ;
reprus= 0.45    ;
repchi= 0.3     ;
repind= 0.75    ;

// ------------------------- endval for nn  -----------------------//

nnadvnam=       0.573599005     ;
nneasnam=       0.110787211     ;
nnmennam=       1.128846267     ;
nnlacnam=       1.091592768     ;
nnjapnam=       0.121486769     ;
nnssanam=       3.212828821     ;
nnrusnam=       0.285373715     ;
nnchinam=       2.348366248     ;
nnindnam=       4.049241021     ;

//-------------------- endval for skilled proportion  -------------------//

phinam= 0.55    ;
phiadv= 0.3     ;
phieas= 0.2     ;
phimen= 0.15    ;
philac= 0.15    ;
phijap= 0.35    ;
phissa= 0.05    ;
phirus= 0.25    ;
phichi= 0.1     ;
phiind= 0.1     ;

//-------------------------- endval for nn  --------------------------//

mmnam= 1 ;
mmadv= 1 ;
mmeas= 1 ;
mmmen= 1 ;
mmlac= 1 ;
mmjap= 1 ;
mmssa= 1 ;
mmrus= 1 ;
mmchi= 1 ;
mmind= 1 ;

// --------------------- endval for consumption tax ------------------//

tcnam=  0.2     ;
tcadv=  0.2     ;
tceas=  0.11    ;
tcmen=  0.11    ;
tclac=  0.11    ;
tcjap=  0.2     ;
tcssa=  0.11    ;
tcrus=  0.11    ;
tcchi=  0.11    ;
tcind=  0.11    ;

//-------------------- endval for country risk rating ----------------//

rankeas=        3.4     ;
rankmen=        3.952380952     ;
ranklac=        5.185185185     ;
rankssa=        6.404761905     ;
rankrus=        6.166666667     ;
rankchi=        3.181818182     ;
rankind=        4.888888889     ;

pimax=0.5;

//-------------------- endval for risk premium pi ---------------------//

pissa=rankssa/7*pimax;
pilac=ranklac/7*pimax;
pimen=rankmen/7*pimax;
pieas=rankeas/7*pimax;
pirus=rankrus/7*pimax;
pichi=rankchi/7*pimax;
piind=rankind/7*pimax;

//----------------------- endval for ratiogdp  ------------------------//

ratiogdpadv= 0.741069062 ;
ratiogdpeas= 0.275195235 ;
ratiogdpmen= 0.166749549 ;
ratiogdplac= 0.21486832 ;
ratiogdpjap= 0.77302212 ;
ratiogdpssa= 0.051913099 ;
ratiogdprus= 0.157734113 ;
ratiogdpchi= 0.128682334 ;
ratiogdpind= 0.074176055 ;

//-------------------------- endval for psi  --------------------------//

psinam= 1       ;
psiadv= 1       ;
psieas= 1       ;
psimen= 1       ;
psilac= 1       ;
psijap= 1       ;
psissa= 1       ;
psirus= 1       ;
psichi= 1       ;
psiind= 1       ;

//-------------------------- endval for tfp  --------------------------//

tfpadv=         0.86190594      ;
tfpchi=         0.12008983      ;
tfpeas=         0.37457905      ;
tfpind=         0.074057072     ;
tfpjap=         0.96311579      ;
tfplac=         0.26885901      ;
tfpmen=         0.1884702       ;
tfprus=         0.24294775      ;
tfpssa=         0.05747975      ;


// -------------------------- non-walras ---------------------------- //

// ------------ endval --------- labs and labu ---------------------- //

labsadv=        1.212431476     ;
labsnam=        2.222791039     ;
labseas=        0.750931459     ;
labsmen=        0.564228589     ;
labslac=        0.553518162     ;
labsjap=        1.412642187     ;
labsssa=        0.170913689     ;
labsrus=        0.901205135     ;
labschi=        0.372999794     ;
labsind=        0.372001999     ;

labuadv=        2.969006776     ;
labunam=        1.908647213     ;
labueas=        3.483725837     ;
labumen=        3.707295336     ;
labulac=        3.646602919     ;
labujap=        2.754465923     ;
labussa=        3.817360082     ;
laburus=        3.153615405     ;
labuchi=        3.896998142     ;
labuind=        3.888017995     ;

// ---------------------- wagu and wags ------------------------- //

wagsadv       =         0.28896;
wagschi       =         0.0669124;
wagseas       =         0.129274;
wagsind       =         0.0380894;
wagsjap       =         0.318184;
wagslac       =         0.107604;
wagsmen       =         0.083527;
wagsnam       =         0.328897;
wagsrus         =       0.064897;
wagsssa         =       0.0285806;

waguadv         =       0.122962;
waguchi         =       0.0223041;
wagueas         =       0.0430914;
waguind         =       0.0126965;
wagujap         =       0.106061;
wagulac         =       0.035868;
wagumen         =       0.0278423;
wagunam         =       0.109632;
wagurus         =       0.0216323;
wagussa         =       0.00952686;

//-------------------- endval for   eta  ------------------------//

etaadv=         0.02798426      ;
etachi=         0.004929291     ;
etaeas=         0.026936141     ;
etaind=         0.011677204     ;
etajap=         0.010524406     ;
etalac= 0.016353353     ;
etamen=         0.025747531     ;
etanam=         0.010414856     ;
etarus=         0.017247165     ;
etassa=         0.036543636     ;


//-------------------- endval for   aa  ------------------------//

aaadv=          0.55367118      ;
aachi=          0.3781675       ;
aaeas=          0.50062749      ;
aaind=          0.37810508      ;
aajap=          0.65077675      ;
aalac=  0.45036167      ;
aamen=          0.43878344      ;
aanam=          0.7699911       ;
aarus=          0.57051699      ;
aassa=          0.27569454      ;


// ---------- initval ------------ wagddu -------------------------------//

wagdduadv=etaadv*wagsadv+(1-etaadv)*waguadv;
wagddunam=etanam*wagsnam+(1-etanam)*wagunam;
wagddussa=etassa*wagsssa+(1-etassa)*wagussa;
wagddulac=etalac*wagslac+(1-etalac)*wagulac;
wagddujap=etajap*wagsjap+(1-etajap)*wagujap;
wagddueas=etaeas*wagseas+(1-etaeas)*wagueas;
wagddumen=etamen*wagsmen+(1-etamen)*wagumen;
wagddurus=etarus*wagsrus+(1-etarus)*wagurus;
wagdduchi=etachi*wagschi+(1-etachi)*waguchi;
wagdduind=etaind*wagsind+(1-etaind)*waguind;

labdduadv       =       2.83422;
labdduchi       =    3.79741;
labddueas       =       3.13659;
labdduind       =       3.66172;
labddujap       =       2.71749;
labddulac       =       3.38703;
labddumen       =       3.33222;
labddunam       =       1.89092;
labddurus       =       2.93369;
labddussa       =       3.10161;

wagddsadv=aaadv*(1-alpha)*kadv^alpha*tfpadv^(1-alpha)*labsadv^(sigma-1)
*(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^((1-alpha-sigma)/sigma);
wagddsnam=aanam*(1-alpha)*knam^alpha*labsnam^(sigma-1)
*(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^((1-alpha-sigma)/sigma);
wagddsssa=aassa*(1-alpha)*kssa^alpha*tfpssa^(1-alpha)*labsssa^(sigma-1)
*(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^((1-alpha-sigma)/sigma);
wagddslac=aalac*(1-alpha)*klac^alpha*tfplac^(1-alpha)*labslac^(sigma-1)
*(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^((1-alpha-sigma)/sigma);
wagddsjap=aajap*(1-alpha)*kjap^alpha*tfpjap^(1-alpha)*labsjap^(sigma-1)
*(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^((1-alpha-sigma)/sigma);
wagddseas=aaeas*(1-alpha)*keas^alpha*tfpeas^(1-alpha)*labseas^(sigma-1)
*(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^((1-alpha-sigma)/sigma);
wagddsmen=aamen*(1-alpha)*kmen^alpha*tfpmen^(1-alpha)*labsmen^(sigma-1)
*(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^((1-alpha-sigma)/sigma);
wagddsrus=aarus*(1-alpha)*krus^alpha*tfprus^(1-alpha)*labsrus^(sigma-1)
*(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^((1-alpha-sigma)/sigma);
wagddschi=aachi*(1-alpha)*kchi^alpha*tfpchi^(1-alpha)*labschi^(sigma-1)
*(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^((1-alpha-sigma)/sigma);
wagddsind=aaind*(1-alpha)*kind^alpha*tfpind^(1-alpha)*labsind^(sigma-1)
*(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^((1-alpha-sigma)/sigma);

// ------------- initval ------- unemeployement ----------------------- //

uradv=(labuadv-labdduadv)/labuadv;
urnam=(labunam-labddunam)/labunam;
urssa=(labussa-labddussa)/labussa;
urlac=(labulac-labddulac)/labulac;
urjap=(labujap-labddujap)/labujap;
ureas=(labueas-labddueas)/labueas;
urmen=(labumen-labddumen)/labumen;
urrus=(laburus-labddurus)/laburus;
urchi=(labuchi-labdduchi)/labuchi;
urind=(labuind-labdduind)/labuind;

// -------------------------------- h ------------------------------- //

hadv=wagsadv/waguadv;
hnam=wagsnam/wagunam;
hssa=wagsssa/wagussa;
hlac=wagslac/wagulac;
hjap=wagsjap/wagujap;
heas=wagseas/wagueas;
hmen=wagsmen/wagumen;
hrus=wagsrus/wagurus;
hchi=wagschi/waguchi;
hind=wagsind/waguind;

// ------------------------- labddu & wagddu ------------------------- //

thetaadv=wagddsadv/wagdduadv;
thetanam=wagddsnam/wagddunam;
thetassa=wagddsssa/wagddussa;
thetalac=wagddslac/wagddulac;
thetajap=wagddsjap/wagddujap;
thetaeas=wagddseas/wagddueas;
thetamen=wagddsmen/wagddumen;
thetarus=wagddsrus/wagddurus;
thetachi=wagddschi/wagdduchi;
thetaind=wagddsind/wagdduind;

// -------------------------------- ces ------------------------------- //

labadv=(aaadv*labsadv^sigma+(1-aaadv)*labdduadv^sigma)^(1/sigma);
labnam=(aanam*labsnam^sigma+(1-aanam)*labddunam^sigma)^(1/sigma);
labssa=(aassa*labsssa^sigma+(1-aassa)*labddussa^sigma)^(1/sigma);
lablac=(aalac*labslac^sigma+(1-aalac)*labddulac^sigma)^(1/sigma);
labjap=(aajap*labsjap^sigma+(1-aajap)*labddujap^sigma)^(1/sigma);
labeas=(aaeas*labseas^sigma+(1-aaeas)*labddueas^sigma)^(1/sigma);
labmen=(aamen*labsmen^sigma+(1-aamen)*labddumen^sigma)^(1/sigma);
labrus=(aarus*labsrus^sigma+(1-aarus)*labddurus^sigma)^(1/sigma);
labchi=(aachi*labschi^sigma+(1-aachi)*labdduchi^sigma)^(1/sigma);
labind=(aaind*labsind^sigma+(1-aaind)*labdduind^sigma)^(1/sigma);

// ---------- end --------- non-walras -------------------------------//

//------------ fin ------------ non-walras ----------------------------//

//--------- new variables --------------------  bengdp ------------//

bengdpadv=(bsadv*((1-lams1adv)*phiadv*(1-edusadv)
+(1-lams2adv)*phiadv*p2adv/mmadv
+(1-lams3adv)*phiadv*p3adv/(mmadv*mmadv)
+(1-lams4adv)*phiadv*p4adv/(mmadv*mmadv*mmadv)
+(1-lams5adv)*phiadv*p5adv/(mmadv*mmadv*mmadv*mmadv)
+(1-lams6adv)*phiadv*p6adv/(mmadv*mmadv*mmadv*mmadv*mmadv)
+(1-lams7adv)*phiadv*p7adv/(mmadv*mmadv*mmadv*mmadv*mmadv*mmadv)
+(1-lams8adv)*phiadv*p8adv/(mmadv*mmadv*mmadv*mmadv*mmadv*mmadv*mmadv))
+buadv*((1-lamu1adv)*(1-phiadv)*(1-eduuadv)
+(1-lamu2adv)*(1-phiadv)*p2adv/mmadv
+(1-lamu3adv)*(1-phiadv)*p3adv/(mmadv*mmadv)
+(1-lamu4adv)*(1-phiadv)*p4adv/(mmadv*mmadv*mmadv)
+(1-lamu5adv)*(1-phiadv)*p5adv/(mmadv*mmadv*mmadv*mmadv)
+(1-lamu6adv)*(1-phiadv)*p6adv/(mmadv*mmadv*mmadv*mmadv*mmadv)
+(1-lamu7adv)*(1-phiadv)*p7adv/(mmadv*mmadv*mmadv*mmadv*mmadv*mmadv)
+(1-lamu8adv)*(1-phiadv)*p8adv/(mmadv*mmadv*mmadv*mmadv*mmadv*mmadv*mmadv)))/gdpadv;

bengdpnam=(bsnam*((1-lams1nam)*phinam*(1-edusnam)
+(1-lams2nam)*phinam*p2nam/mmnam
+(1-lams3nam)*phinam*p3nam/(mmnam*mmnam)
+(1-lams4nam)*phinam*p4nam/(mmnam*mmnam*mmnam)
+(1-lams5nam)*phinam*p5nam/(mmnam*mmnam*mmnam*mmnam)
+(1-lams6nam)*phinam*p6nam/(mmnam*mmnam*mmnam*mmnam*mmnam)
+(1-lams7nam)*phinam*p7nam/(mmnam*mmnam*mmnam*mmnam*mmnam*mmnam)
+(1-lams8nam)*phinam*p8nam/(mmnam*mmnam*mmnam*mmnam*mmnam*mmnam*mmnam))
+bunam*((1-lamu1nam)*(1-phinam)*(1-eduunam)
+(1-lamu2nam)*(1-phinam)*p2nam/mmnam
+(1-lamu3nam)*(1-phinam)*p3nam/(mmnam*mmnam)
+(1-lamu4nam)*(1-phinam)*p4nam/(mmnam*mmnam*mmnam)
+(1-lamu5nam)*(1-phinam)*p5nam/(mmnam*mmnam*mmnam*mmnam)
+(1-lamu6nam)*(1-phinam)*p6nam/(mmnam*mmnam*mmnam*mmnam*mmnam)
+(1-lamu7nam)*(1-phinam)*p7nam/(mmnam*mmnam*mmnam*mmnam*mmnam*mmnam)
+(1-lamu8nam)*(1-phinam)*p8nam/(mmnam*mmnam*mmnam*mmnam*mmnam*mmnam*mmnam)))/gdpnam;

bengdpssa=(bsssa*((1-lams1ssa)*phissa*(1-edusssa)
+(1-lams2ssa)*phissa*p2ssa/mmssa
+(1-lams3ssa)*phissa*p3ssa/(mmssa*mmssa)
+(1-lams4ssa)*phissa*p4ssa/(mmssa*mmssa*mmssa)
+(1-lams5ssa)*phissa*p5ssa/(mmssa*mmssa*mmssa*mmssa)
+(1-lams6ssa)*phissa*p6ssa/(mmssa*mmssa*mmssa*mmssa*mmssa)
+(1-lams7ssa)*phissa*p7ssa/(mmssa*mmssa*mmssa*mmssa*mmssa*mmssa)
+(1-lams8ssa)*phissa*p8ssa/(mmssa*mmssa*mmssa*mmssa*mmssa*mmssa*mmssa))
+bussa*((1-lamu1ssa)*(1-phissa)*(1-eduussa)
+(1-lamu2ssa)*(1-phissa)*p2ssa/mmssa
+(1-lamu3ssa)*(1-phissa)*p3ssa/(mmssa*mmssa)
+(1-lamu4ssa)*(1-phissa)*p4ssa/(mmssa*mmssa*mmssa)
+(1-lamu5ssa)*(1-phissa)*p5ssa/(mmssa*mmssa*mmssa*mmssa)
+(1-lamu6ssa)*(1-phissa)*p6ssa/(mmssa*mmssa*mmssa*mmssa*mmssa)
+(1-lamu7ssa)*(1-phissa)*p7ssa/(mmssa*mmssa*mmssa*mmssa*mmssa*mmssa)
+(1-lamu8ssa)*(1-phissa)*p8ssa/(mmssa*mmssa*mmssa*mmssa*mmssa*mmssa*mmssa)))/gdpssa;

bengdplac=(bslac*((1-lams1lac)*philac*(1-eduslac)
+(1-lams2lac)*philac*p2lac/mmlac
+(1-lams3lac)*philac*p3lac/(mmlac*mmlac)
+(1-lams4lac)*philac*p4lac/(mmlac*mmlac*mmlac)
+(1-lams5lac)*philac*p5lac/(mmlac*mmlac*mmlac*mmlac)
+(1-lams6lac)*philac*p6lac/(mmlac*mmlac*mmlac*mmlac*mmlac)
+(1-lams7lac)*philac*p7lac/(mmlac*mmlac*mmlac*mmlac*mmlac*mmlac)
+(1-lams8lac)*philac*p8lac/(mmlac*mmlac*mmlac*mmlac*mmlac*mmlac*mmlac))
+bulac*((1-lamu1lac)*(1-philac)*(1-eduulac)
+(1-lamu2lac)*(1-philac)*p2lac/mmlac
+(1-lamu3lac)*(1-philac)*p3lac/(mmlac*mmlac)
+(1-lamu4lac)*(1-philac)*p4lac/(mmlac*mmlac*mmlac)
+(1-lamu5lac)*(1-philac)*p5lac/(mmlac*mmlac*mmlac*mmlac)
+(1-lamu6lac)*(1-philac)*p6lac/(mmlac*mmlac*mmlac*mmlac*mmlac)
+(1-lamu7lac)*(1-philac)*p7lac/(mmlac*mmlac*mmlac*mmlac*mmlac*mmlac)
+(1-lamu8lac)*(1-philac)*p8lac/(mmlac*mmlac*mmlac*mmlac*mmlac*mmlac*mmlac)))/gdplac;

bengdpjap=(bsjap*((1-lams1jap)*phijap*(1-edusjap)
+(1-lams2jap)*phijap*p2jap/mmjap
+(1-lams3jap)*phijap*p3jap/(mmjap*mmjap)
+(1-lams4jap)*phijap*p4jap/(mmjap*mmjap*mmjap)
+(1-lams5jap)*phijap*p5jap/(mmjap*mmjap*mmjap*mmjap)
+(1-lams6jap)*phijap*p6jap/(mmjap*mmjap*mmjap*mmjap*mmjap)
+(1-lams7jap)*phijap*p7jap/(mmjap*mmjap*mmjap*mmjap*mmjap*mmjap)
+(1-lams8jap)*phijap*p8jap/(mmjap*mmjap*mmjap*mmjap*mmjap*mmjap*mmjap))
+bujap*((1-lamu1jap)*(1-phijap)*(1-eduujap)
+(1-lamu2jap)*(1-phijap)*p2jap/mmjap
+(1-lamu3jap)*(1-phijap)*p3jap/(mmjap*mmjap)
+(1-lamu4jap)*(1-phijap)*p4jap/(mmjap*mmjap*mmjap)
+(1-lamu5jap)*(1-phijap)*p5jap/(mmjap*mmjap*mmjap*mmjap)
+(1-lamu6jap)*(1-phijap)*p6jap/(mmjap*mmjap*mmjap*mmjap*mmjap)
+(1-lamu7jap)*(1-phijap)*p7jap/(mmjap*mmjap*mmjap*mmjap*mmjap*mmjap)
+(1-lamu8jap)*(1-phijap)*p8jap/(mmjap*mmjap*mmjap*mmjap*mmjap*mmjap*mmjap)))/gdpjap;

bengdprus=(bsrus*((1-lams1rus)*phirus*(1-edusrus)
+(1-lams2rus)*phirus*p2rus/mmrus
+(1-lams3rus)*phirus*p3rus/(mmrus*mmrus)
+(1-lams4rus)*phirus*p4rus/(mmrus*mmrus*mmrus)
+(1-lams5rus)*phirus*p5rus/(mmrus*mmrus*mmrus*mmrus)
+(1-lams6rus)*phirus*p6rus/(mmrus*mmrus*mmrus*mmrus*mmrus)
+(1-lams7rus)*phirus*p7rus/(mmrus*mmrus*mmrus*mmrus*mmrus*mmrus)
+(1-lams8rus)*phirus*p8rus/(mmrus*mmrus*mmrus*mmrus*mmrus*mmrus*mmrus))
+burus*((1-lamu1rus)*(1-phirus)*(1-eduurus)
+(1-lamu2rus)*(1-phirus)*p2rus/mmrus
+(1-lamu3rus)*(1-phirus)*p3rus/(mmrus*mmrus)
+(1-lamu4rus)*(1-phirus)*p4rus/(mmrus*mmrus*mmrus)
+(1-lamu5rus)*(1-phirus)*p5rus/(mmrus*mmrus*mmrus*mmrus)
+(1-lamu6rus)*(1-phirus)*p6rus/(mmrus*mmrus*mmrus*mmrus*mmrus)
+(1-lamu7rus)*(1-phirus)*p7rus/(mmrus*mmrus*mmrus*mmrus*mmrus*mmrus)
+(1-lamu8rus)*(1-phirus)*p8rus/(mmrus*mmrus*mmrus*mmrus*mmrus*mmrus*mmrus)))/gdprus;

bengdpmen=(bsmen*((1-lams1men)*phimen*(1-edusmen)
+(1-lams2men)*phimen*p2men/mmmen
+(1-lams3men)*phimen*p3men/(mmmen*mmmen)
+(1-lams4men)*phimen*p4men/(mmmen*mmmen*mmmen)
+(1-lams5men)*phimen*p5men/(mmmen*mmmen*mmmen*mmmen)
+(1-lams6men)*phimen*p6men/(mmmen*mmmen*mmmen*mmmen*mmmen)
+(1-lams7men)*phimen*p7men/(mmmen*mmmen*mmmen*mmmen*mmmen*mmmen)
+(1-lams8men)*phimen*p8men/(mmmen*mmmen*mmmen*mmmen*mmmen*mmmen*mmmen))
+bumen*((1-lamu1men)*(1-phimen)*(1-eduumen)
+(1-lamu2men)*(1-phimen)*p2men/mmmen
+(1-lamu3men)*(1-phimen)*p3men/(mmmen*mmmen)
+(1-lamu4men)*(1-phimen)*p4men/(mmmen*mmmen*mmmen)
+(1-lamu5men)*(1-phimen)*p5men/(mmmen*mmmen*mmmen*mmmen)
+(1-lamu6men)*(1-phimen)*p6men/(mmmen*mmmen*mmmen*mmmen*mmmen)
+(1-lamu7men)*(1-phimen)*p7men/(mmmen*mmmen*mmmen*mmmen*mmmen*mmmen)
+(1-lamu8men)*(1-phimen)*p8men/(mmmen*mmmen*mmmen*mmmen*mmmen*mmmen*mmmen)))/gdpmen;

bengdpeas=(bseas*((1-lams1eas)*phieas*(1-eduseas)
+(1-lams2eas)*phieas*p2eas/mmeas
+(1-lams3eas)*phieas*p3eas/(mmeas*mmeas)
+(1-lams4eas)*phieas*p4eas/(mmeas*mmeas*mmeas)
+(1-lams5eas)*phieas*p5eas/(mmeas*mmeas*mmeas*mmeas)
+(1-lams6eas)*phieas*p6eas/(mmeas*mmeas*mmeas*mmeas*mmeas)
+(1-lams7eas)*phieas*p7eas/(mmeas*mmeas*mmeas*mmeas*mmeas*mmeas)
+(1-lams8eas)*phieas*p8eas/(mmeas*mmeas*mmeas*mmeas*mmeas*mmeas*mmeas))
+bueas*((1-lamu1eas)*(1-phieas)*(1-eduueas)
+(1-lamu2eas)*(1-phieas)*p2eas/mmeas
+(1-lamu3eas)*(1-phieas)*p3eas/(mmeas*mmeas)
+(1-lamu4eas)*(1-phieas)*p4eas/(mmeas*mmeas*mmeas)
+(1-lamu5eas)*(1-phieas)*p5eas/(mmeas*mmeas*mmeas*mmeas)
+(1-lamu6eas)*(1-phieas)*p6eas/(mmeas*mmeas*mmeas*mmeas*mmeas)
+(1-lamu7eas)*(1-phieas)*p7eas/(mmeas*mmeas*mmeas*mmeas*mmeas*mmeas)
+(1-lamu8eas)*(1-phieas)*p8eas/(mmeas*mmeas*mmeas*mmeas*mmeas*mmeas*mmeas)))/gdpeas;

bengdpchi=(bschi*((1-lams1chi)*phichi*(1-eduschi)
+(1-lams2chi)*phichi*p2chi/mmchi
+(1-lams3chi)*phichi*p3chi/(mmchi*mmchi)
+(1-lams4chi)*phichi*p4chi/(mmchi*mmchi*mmchi)
+(1-lams5chi)*phichi*p5chi/(mmchi*mmchi*mmchi*mmchi)
+(1-lams6chi)*phichi*p6chi/(mmchi*mmchi*mmchi*mmchi*mmchi)
+(1-lams7chi)*phichi*p7chi/(mmchi*mmchi*mmchi*mmchi*mmchi*mmchi)
+(1-lams8chi)*phichi*p8chi/(mmchi*mmchi*mmchi*mmchi*mmchi*mmchi*mmchi))
+buchi*((1-lamu1chi)*(1-phichi)*(1-eduuchi)
+(1-lamu2chi)*(1-phichi)*p2chi/mmchi
+(1-lamu3chi)*(1-phichi)*p3chi/(mmchi*mmchi)
+(1-lamu4chi)*(1-phichi)*p4chi/(mmchi*mmchi*mmchi)
+(1-lamu5chi)*(1-phichi)*p5chi/(mmchi*mmchi*mmchi*mmchi)
+(1-lamu6chi)*(1-phichi)*p6chi/(mmchi*mmchi*mmchi*mmchi*mmchi)
+(1-lamu7chi)*(1-phichi)*p7chi/(mmchi*mmchi*mmchi*mmchi*mmchi*mmchi)
+(1-lamu8chi)*(1-phichi)*p8chi/(mmchi*mmchi*mmchi*mmchi*mmchi*mmchi*mmchi)))/gdpchi;

bengdpind=(bsind*((1-lams1ind)*phiind*(1-edusind)
+(1-lams2ind)*phiind*p2ind/mmind
+(1-lams3ind)*phiind*p3ind/(mmind*mmind)
+(1-lams4ind)*phiind*p4ind/(mmind*mmind*mmind)
+(1-lams5ind)*phiind*p5ind/(mmind*mmind*mmind*mmind)
+(1-lams6ind)*phiind*p6ind/(mmind*mmind*mmind*mmind*mmind)
+(1-lams7ind)*phiind*p7ind/(mmind*mmind*mmind*mmind*mmind*mmind)
+(1-lams8ind)*phiind*p8ind/(mmind*mmind*mmind*mmind*mmind*mmind*mmind))
+buind*((1-lamu1ind)*(1-phiind)*(1-eduuind)
+(1-lamu2ind)*(1-phiind)*p2ind/mmind
+(1-lamu3ind)*(1-phiind)*p3ind/(mmind*mmind)
+(1-lamu4ind)*(1-phiind)*p4ind/(mmind*mmind*mmind)
+(1-lamu5ind)*(1-phiind)*p5ind/(mmind*mmind*mmind*mmind)
+(1-lamu6ind)*(1-phiind)*p6ind/(mmind*mmind*mmind*mmind*mmind)
+(1-lamu7ind)*(1-phiind)*p7ind/(mmind*mmind*mmind*mmind*mmind*mmind)
+(1-lamu8ind)*(1-phiind)*p8ind/(mmind*mmind*mmind*mmind*mmind*mmind*mmind)))/gdpind;

// --------------- endval --------- sup ------------------------------- //

supadv     =            2.33651;
supchi      =           3.18865;
supeas       =          3.26033;
supind        =         3.63779;
supjap        =         2.25815;
suplac        =         3.07692;
supmen        =         3.32751;
supnam        =         2.38249;
suprus        =         4.0982;
supssa        =         5.66508;

// ------------------ endval ------ trgdp ----------------------------- //

trgdpadv=(psiadv*wagddsadv*(tr1sadv*phiadv
+tr2sadv*phiadv*p2adv/mmadv
+tr3sadv*phiadv*p3adv/(mmadv^2)
+tr4sadv*phiadv*p4adv/(mmadv^3)
+tr5sadv*phiadv*p5adv/(mmadv^4)
+tr6sadv*phiadv*p6adv/(mmadv^5)
+tr7sadv*phiadv*p7adv/(mmadv^6)
+tr8sadv*phiadv*p8adv/(mmadv^7))
+psiadv*waguadv*(tr1uadv*(1-phiadv)
+tr2uadv*(1-phiadv)*p2adv/mmadv
+tr3uadv*(1-phiadv)*p3adv/(mmadv^2)
+tr4uadv*(1-phiadv)*p4adv/(mmadv^3)
+tr5uadv*(1-phiadv)*p5adv/(mmadv^4)
+tr6uadv*(1-phiadv)*p6adv/(mmadv^5)
+tr7uadv*(1-phiadv)*p7adv/(mmadv^6)
+tr8uadv*(1-phiadv)*p8adv/(mmadv^7)))/gdpadv;

trgdpnam=(psinam*wagddsnam*(tr1snam*phinam
+tr2snam*phinam*p2nam/mmnam
+tr3snam*phinam*p3nam/(mmnam^2)
+tr4snam*phinam*p4nam/(mmnam^3)
+tr5snam*phinam*p5nam/(mmnam^4)
+tr6snam*phinam*p6nam/(mmnam^5)
+tr7snam*phinam*p7nam/(mmnam^6)
+tr8snam*phinam*p8nam/(mmnam^7))
+psinam*wagunam*(tr1unam*(1-phinam)
+tr2unam*(1-phinam)*p2nam/mmnam
+tr3unam*(1-phinam)*p3nam/(mmnam^2)
+tr4unam*(1-phinam)*p4nam/(mmnam^3)
+tr5unam*(1-phinam)*p5nam/(mmnam^4)
+tr6unam*(1-phinam)*p6nam/(mmnam^5)
+tr7unam*(1-phinam)*p7nam/(mmnam^6)
+tr8unam*(1-phinam)*p8nam/(mmnam^7)))/gdpnam;

trgdpssa=(psissa*wagddsssa*(tr1sssa*phissa
+tr2sssa*phissa*p2ssa/mmssa
+tr3sssa*phissa*p3ssa/(mmssa^2)
+tr4sssa*phissa*p4ssa/(mmssa^3)
+tr5sssa*phissa*p5ssa/(mmssa^4)
+tr6sssa*phissa*p6ssa/(mmssa^5)
+tr7sssa*phissa*p7ssa/(mmssa^6)
+tr8sssa*phissa*p8ssa/(mmssa^7))
+psissa*wagussa*(tr1ussa*(1-phissa)
+tr2ussa*(1-phissa)*p2ssa/mmssa
+tr3ussa*(1-phissa)*p3ssa/(mmssa^2)
+tr4ussa*(1-phissa)*p4ssa/(mmssa^3)
+tr5ussa*(1-phissa)*p5ssa/(mmssa^4)
+tr6ussa*(1-phissa)*p6ssa/(mmssa^5)
+tr7ussa*(1-phissa)*p7ssa/(mmssa^6)
+tr8ussa*(1-phissa)*p8ssa/(mmssa^7)))/gdpssa;

trgdplac=(psilac*wagddslac*(tr1slac*philac
+tr2slac*philac*p2lac/mmlac
+tr3slac*philac*p3lac/(mmlac^2)
+tr4slac*philac*p4lac/(mmlac^3)
+tr5slac*philac*p5lac/(mmlac^4)
+tr6slac*philac*p6lac/(mmlac^5)
+tr7slac*philac*p7lac/(mmlac^6)
+tr8slac*philac*p8lac/(mmlac^7))
+psilac*wagulac*(tr1ulac*(1-philac)
+tr2ulac*(1-philac)*p2lac/mmlac
+tr3ulac*(1-philac)*p3lac/(mmlac^2)
+tr4ulac*(1-philac)*p4lac/(mmlac^3)
+tr5ulac*(1-philac)*p5lac/(mmlac^4)
+tr6ulac*(1-philac)*p6lac/(mmlac^5)
+tr7ulac*(1-philac)*p7lac/(mmlac^6)
+tr8ulac*(1-philac)*p8lac/(mmlac^7)))/gdplac;

trgdpjap=(psijap*wagddsjap*(tr1sjap*phijap
+tr2sjap*phijap*p2jap/mmjap
+tr3sjap*phijap*p3jap/(mmjap^2)
+tr4sjap*phijap*p4jap/(mmjap^3)
+tr5sjap*phijap*p5jap/(mmjap^4)
+tr6sjap*phijap*p6jap/(mmjap^5)
+tr7sjap*phijap*p7jap/(mmjap^6)
+tr8sjap*phijap*p8jap/(mmjap^7))
+psijap*wagujap*(tr1ujap*(1-phijap)
+tr2ujap*(1-phijap)*p2jap/mmjap
+tr3ujap*(1-phijap)*p3jap/(mmjap^2)
+tr4ujap*(1-phijap)*p4jap/(mmjap^3)
+tr5ujap*(1-phijap)*p5jap/(mmjap^4)
+tr6ujap*(1-phijap)*p6jap/(mmjap^5)
+tr7ujap*(1-phijap)*p7jap/(mmjap^6)
+tr8ujap*(1-phijap)*p8jap/(mmjap^7)))/gdpjap;

trgdprus=(psirus*wagddsrus*(tr1srus*phirus
+tr2srus*phirus*p2rus/mmrus
+tr3srus*phirus*p3rus/(mmrus^2)
+tr4srus*phirus*p4rus/(mmrus^3)
+tr5srus*phirus*p5rus/(mmrus^4)
+tr6srus*phirus*p6rus/(mmrus^5)
+tr7srus*phirus*p7rus/(mmrus^6)
+tr8srus*phirus*p8rus/(mmrus^7))
+psirus*wagurus*(tr1urus*(1-phirus)
+tr2urus*(1-phirus)*p2rus/mmrus
+tr3urus*(1-phirus)*p3rus/(mmrus^2)
+tr4urus*(1-phirus)*p4rus/(mmrus^3)
+tr5urus*(1-phirus)*p5rus/(mmrus^4)
+tr6urus*(1-phirus)*p6rus/(mmrus^5)
+tr7urus*(1-phirus)*p7rus/(mmrus^6)
+tr8urus*(1-phirus)*p8rus/(mmrus^7)))/gdprus;

trgdpmen=(psimen*wagddsmen*(tr1smen*phimen
+tr2smen*phimen*p2men/mmmen
+tr3smen*phimen*p3men/(mmmen^2)
+tr4smen*phimen*p4men/(mmmen^3)
+tr5smen*phimen*p5men/(mmmen^4)
+tr6smen*phimen*p6men/(mmmen^5)
+tr7smen*phimen*p7men/(mmmen^6)
+tr8smen*phimen*p8men/(mmmen^7))
+psimen*wagumen*(tr1umen*(1-phimen)
+tr2umen*(1-phimen)*p2men/mmmen
+tr3umen*(1-phimen)*p3men/(mmmen^2)
+tr4umen*(1-phimen)*p4men/(mmmen^3)
+tr5umen*(1-phimen)*p5men/(mmmen^4)
+tr6umen*(1-phimen)*p6men/(mmmen^5)
+tr7umen*(1-phimen)*p7men/(mmmen^6)
+tr8umen*(1-phimen)*p8men/(mmmen^7)))/gdpmen;

trgdpeas=(psieas*wagddseas*(tr1seas*phieas
+tr2seas*phieas*p2eas/mmeas
+tr3seas*phieas*p3eas/(mmeas^2)
+tr4seas*phieas*p4eas/(mmeas^3)
+tr5seas*phieas*p5eas/(mmeas^4)
+tr6seas*phieas*p6eas/(mmeas^5)
+tr7seas*phieas*p7eas/(mmeas^6)
+tr8seas*phieas*p8eas/(mmeas^7))
+psieas*wagueas*(tr1ueas*(1-phieas)
+tr2ueas*(1-phieas)*p2eas/mmeas
+tr3ueas*(1-phieas)*p3eas/(mmeas^2)
+tr4ueas*(1-phieas)*p4eas/(mmeas^3)
+tr5ueas*(1-phieas)*p5eas/(mmeas^4)
+tr6ueas*(1-phieas)*p6eas/(mmeas^5)
+tr7ueas*(1-phieas)*p7eas/(mmeas^6)
+tr8ueas*(1-phieas)*p8eas/(mmeas^7)))/gdpeas;

trgdpchi=(psichi*wagddschi*(tr1schi*phichi
+tr2schi*phichi*p2chi/mmchi
+tr3schi*phichi*p3chi/(mmchi^2)
+tr4schi*phichi*p4chi/(mmchi^3)
+tr5schi*phichi*p5chi/(mmchi^4)
+tr6schi*phichi*p6chi/(mmchi^5)
+tr7schi*phichi*p7chi/(mmchi^6)
+tr8schi*phichi*p8chi/(mmchi^7))
+psichi*waguchi*(tr1uchi*(1-phichi)
+tr2uchi*(1-phichi)*p2chi/mmchi
+tr3uchi*(1-phichi)*p3chi/(mmchi^2)
+tr4uchi*(1-phichi)*p4chi/(mmchi^3)
+tr5uchi*(1-phichi)*p5chi/(mmchi^4)
+tr6uchi*(1-phichi)*p6chi/(mmchi^5)
+tr7uchi*(1-phichi)*p7chi/(mmchi^6)
+tr8uchi*(1-phichi)*p8chi/(mmchi^7)))/gdpchi;

trgdpind=(psiind*wagddsind*(tr1sind*phiind
+tr2sind*phiind*p2ind/mmind
+tr3sind*phiind*p3ind/(mmind^2)
+tr4sind*phiind*p4ind/(mmind^3)
+tr5sind*phiind*p5ind/(mmind^4)
+tr6sind*phiind*p6ind/(mmind^5)
+tr7sind*phiind*p7ind/(mmind^6)
+tr8sind*phiind*p8ind/(mmind^7))
+psiind*waguind*(tr1uind*(1-phiind)
+tr2uind*(1-phiind)*p2ind/mmind
+tr3uind*(1-phiind)*p3ind/(mmind^2)
+tr4uind*(1-phiind)*p4ind/(mmind^3)
+tr5uind*(1-phiind)*p5ind/(mmind^4)
+tr6uind*(1-phiind)*p6ind/(mmind^5)
+tr7uind*(1-phiind)*p7ind/(mmind^6)
+tr8uind*(1-phiind)*p8ind/(mmind^7)))/gdpind;

// --------- endval ---------- rapport lab-lab sans h --------------- //

lablabsshadv=labadv/((1-edusadv)*phiadv*lams1adv
+phiadv*lams2adv*p2adv/mmadv
+phiadv*lams3adv*p3adv/(mmadv^2)
+phiadv*lams4adv*p4adv/(mmadv^3)
+phiadv*lams5adv*p5adv/(mmadv^4)
+phiadv*lams6adv*p6adv/(mmadv^5)
+phiadv*lams7adv*p7adv/(mmadv^6)
+phiadv*lams8adv*p8adv/(mmadv^7)
+(1-eduuadv)*(1-phiadv)*lamu1adv
+(1-phiadv)*lamu2adv*p2adv/mmadv
+(1-phiadv)*lamu3adv*p3adv/(mmadv^2)
+(1-phiadv)*lamu4adv*p4adv/(mmadv^3)
+(1-phiadv)*lamu5adv*p5adv/(mmadv^4)
+(1-phiadv)*lamu6adv*p6adv/(mmadv^5)
+(1-phiadv)*lamu7adv*p7adv/(mmadv^6)
+(1-phiadv)*lamu8adv*p8adv/(mmadv^7));

lablabsshnam=labnam/((1-edusnam)*phinam*lams1nam
+phinam*lams2nam*p2nam/mmnam
+phinam*lams3nam*p3nam/(mmnam^2)
+phinam*lams4nam*p4nam/(mmnam^3)
+phinam*lams5nam*p5nam/(mmnam^4)
+phinam*lams6nam*p6nam/(mmnam^5)
+phinam*lams7nam*p7nam/(mmnam^6)
+phinam*lams8nam*p8nam/(mmnam^7)
+(1-eduunam)*(1-phinam)*lamu1nam
+(1-phinam)*lamu2nam*p2nam/mmnam
+(1-phinam)*lamu3nam*p3nam/(mmnam^2)
+(1-phinam)*lamu4nam*p4nam/(mmnam^3)
+(1-phinam)*lamu5nam*p5nam/(mmnam^4)
+(1-phinam)*lamu6nam*p6nam/(mmnam^5)
+(1-phinam)*lamu7nam*p7nam/(mmnam^6)
+(1-phinam)*lamu8nam*p8nam/(mmnam^7));

lablabsshjap=labjap/((1-edusjap)*phijap*lams1jap
+phijap*lams2jap*p2jap/mmjap
+phijap*lams3jap*p3jap/(mmjap^2)
+phijap*lams4jap*p4jap/(mmjap^3)
+phijap*lams5jap*p5jap/(mmjap^4)
+phijap*lams6jap*p6jap/(mmjap^5)
+phijap*lams7jap*p7jap/(mmjap^6)
+phijap*lams8jap*p8jap/(mmjap^7)
+(1-eduujap)*(1-phijap)*lamu1jap
+(1-phijap)*lamu2jap*p2jap/mmjap
+(1-phijap)*lamu3jap*p3jap/(mmjap^2)
+(1-phijap)*lamu4jap*p4jap/(mmjap^3)
+(1-phijap)*lamu5jap*p5jap/(mmjap^4)
+(1-phijap)*lamu6jap*p6jap/(mmjap^5)
+(1-phijap)*lamu7jap*p7jap/(mmjap^6)
+(1-phijap)*lamu8jap*p8jap/(mmjap^7));

lablabsshlac=lablac/((1-eduslac)*philac*lams1lac
+philac*lams2lac*p2lac/mmlac
+philac*lams3lac*p3lac/(mmlac^2)
+philac*lams4lac*p4lac/(mmlac^3)
+philac*lams5lac*p5lac/(mmlac^4)
+philac*lams6lac*p6lac/(mmlac^5)
+philac*lams7lac*p7lac/(mmlac^6)
+philac*lams8lac*p8lac/(mmlac^7)
+(1-eduulac)*(1-philac)*lamu1lac
+(1-philac)*lamu2lac*p2lac/mmlac
+(1-philac)*lamu3lac*p3lac/(mmlac^2)
+(1-philac)*lamu4lac*p4lac/(mmlac^3)
+(1-philac)*lamu5lac*p5lac/(mmlac^4)
+(1-philac)*lamu6lac*p6lac/(mmlac^5)
+(1-philac)*lamu7lac*p7lac/(mmlac^6)
+(1-philac)*lamu8lac*p8lac/(mmlac^7));

lablabssheas=labeas/((1-eduseas)*phieas*lams1eas
+phieas*lams2eas*p2eas/mmeas
+phieas*lams3eas*p3eas/(mmeas^2)
+phieas*lams4eas*p4eas/(mmeas^3)
+phieas*lams5eas*p5eas/(mmeas^4)
+phieas*lams6eas*p6eas/(mmeas^5)
+phieas*lams7eas*p7eas/(mmeas^6)
+phieas*lams8eas*p8eas/(mmeas^7)
+(1-eduueas)*(1-phieas)*lamu1eas
+(1-phieas)*lamu2eas*p2eas/mmeas
+(1-phieas)*lamu3eas*p3eas/(mmeas^2)
+(1-phieas)*lamu4eas*p4eas/(mmeas^3)
+(1-phieas)*lamu5eas*p5eas/(mmeas^4)
+(1-phieas)*lamu6eas*p6eas/(mmeas^5)
+(1-phieas)*lamu7eas*p7eas/(mmeas^6)
+(1-phieas)*lamu8eas*p8eas/(mmeas^7));

lablabsshmen=labmen/((1-edusmen)*phimen*lams1men
+phimen*lams2men*p2men/mmmen
+phimen*lams3men*p3men/(mmmen^2)
+phimen*lams4men*p4men/(mmmen^3)
+phimen*lams5men*p5men/(mmmen^4)
+phimen*lams6men*p6men/(mmmen^5)
+phimen*lams7men*p7men/(mmmen^6)
+phimen*lams8men*p8men/(mmmen^7)
+(1-eduumen)*(1-phimen)*lamu1men
+(1-phimen)*lamu2men*p2men/mmmen
+(1-phimen)*lamu3men*p3men/(mmmen^2)
+(1-phimen)*lamu4men*p4men/(mmmen^3)
+(1-phimen)*lamu5men*p5men/(mmmen^4)
+(1-phimen)*lamu6men*p6men/(mmmen^5)
+(1-phimen)*lamu7men*p7men/(mmmen^6)
+(1-phimen)*lamu8men*p8men/(mmmen^7));

lablabsshrus=labrus/((1-edusrus)*phirus*lams1rus
+phirus*lams2rus*p2rus/mmrus
+phirus*lams3rus*p3rus/(mmrus^2)
+phirus*lams4rus*p4rus/(mmrus^3)
+phirus*lams5rus*p5rus/(mmrus^4)
+phirus*lams6rus*p6rus/(mmrus^5)
+phirus*lams7rus*p7rus/(mmrus^6)
+phirus*lams8rus*p8rus/(mmrus^7)
+(1-eduurus)*(1-phirus)*lamu1rus
+(1-phirus)*lamu2rus*p2rus/mmrus
+(1-phirus)*lamu3rus*p3rus/(mmrus^2)
+(1-phirus)*lamu4rus*p4rus/(mmrus^3)
+(1-phirus)*lamu5rus*p5rus/(mmrus^4)
+(1-phirus)*lamu6rus*p6rus/(mmrus^5)
+(1-phirus)*lamu7rus*p7rus/(mmrus^6)
+(1-phirus)*lamu8rus*p8rus/(mmrus^7));

lablabsshind=labind/((1-edusind)*phiind*lams1ind
+phiind*lams2ind*p2ind/mmind
+phiind*lams3ind*p3ind/(mmind^2)
+phiind*lams4ind*p4ind/(mmind^3)
+phiind*lams5ind*p5ind/(mmind^4)
+phiind*lams6ind*p6ind/(mmind^5)
+phiind*lams7ind*p7ind/(mmind^6)
+phiind*lams8ind*p8ind/(mmind^7)
+(1-eduuind)*(1-phiind)*lamu1ind
+(1-phiind)*lamu2ind*p2ind/mmind
+(1-phiind)*lamu3ind*p3ind/(mmind^2)
+(1-phiind)*lamu4ind*p4ind/(mmind^3)
+(1-phiind)*lamu5ind*p5ind/(mmind^4)
+(1-phiind)*lamu6ind*p6ind/(mmind^5)
+(1-phiind)*lamu7ind*p7ind/(mmind^6)
+(1-phiind)*lamu8ind*p8ind/(mmind^7));

lablabsshchi=labchi/((1-eduschi)*phichi*lams1chi
+phichi*lams2chi*p2chi/mmchi
+phichi*lams3chi*p3chi/(mmchi^2)
+phichi*lams4chi*p4chi/(mmchi^3)
+phichi*lams5chi*p5chi/(mmchi^4)
+phichi*lams6chi*p6chi/(mmchi^5)
+phichi*lams7chi*p7chi/(mmchi^6)
+phichi*lams8chi*p8chi/(mmchi^7)
+(1-eduuchi)*(1-phichi)*lamu1chi
+(1-phichi)*lamu2chi*p2chi/mmchi
+(1-phichi)*lamu3chi*p3chi/(mmchi^2)
+(1-phichi)*lamu4chi*p4chi/(mmchi^3)
+(1-phichi)*lamu5chi*p5chi/(mmchi^4)
+(1-phichi)*lamu6chi*p6chi/(mmchi^5)
+(1-phichi)*lamu7chi*p7chi/(mmchi^6)
+(1-phichi)*lamu8chi*p8chi/(mmchi^7));

lablabsshssa=labssa/((1-edusssa)*phissa*lams1ssa
+phissa*lams2ssa*p2ssa/mmssa
+phissa*lams3ssa*p3ssa/(mmssa^2)
+phissa*lams4ssa*p4ssa/(mmssa^3)
+phissa*lams5ssa*p5ssa/(mmssa^4)
+phissa*lams6ssa*p6ssa/(mmssa^5)
+phissa*lams7ssa*p7ssa/(mmssa^6)
+phissa*lams8ssa*p8ssa/(mmssa^7)
+(1-eduussa)*(1-phissa)*lamu1ssa
+(1-phissa)*lamu2ssa*p2ssa/mmssa
+(1-phissa)*lamu3ssa*p3ssa/(mmssa^2)
+(1-phissa)*lamu4ssa*p4ssa/(mmssa^3)
+(1-phissa)*lamu5ssa*p5ssa/(mmssa^4)
+(1-phissa)*lamu6ssa*p6ssa/(mmssa^5)
+(1-phissa)*lamu7ssa*p7ssa/(mmssa^6)
+(1-phissa)*lamu8ssa*p8ssa/(mmssa^7));

// --------- endval ---------- tx de prélèvement --------------- //

taxesgdpadv=(tauadv*(wagddsadv*labsadv+wagdduadv*labdduadv)
+tcadv*(phiadv*xs1adv
+xs2adv*phiadv*p2adv/mmadv
+xs3adv*phiadv*p3adv/(mmadv^2)
+xs4adv*phiadv*p4adv/(mmadv^3)
+xs5adv*phiadv*p5adv/(mmadv^4)
+xs6adv*phiadv*p6adv/(mmadv^5)
+xs7adv*phiadv*p7adv/(mmadv^6)
+xs8adv*phiadv*p8adv/(mmadv^7))
+tcadv*((1-phiadv)*xu1adv
+xu2adv*(1-phiadv)*p2adv/mmadv
+xu3adv*(1-phiadv)*p3adv/(mmadv^2)
+xu4adv*(1-phiadv)*p4adv/(mmadv^3)
+xu5adv*(1-phiadv)*p5adv/(mmadv^4)
+xu6adv*(1-phiadv)*p6adv/(mmadv^5)
+xu7adv*(1-phiadv)*p7adv/(mmadv^6)
+xu8adv*(1-phiadv)*p8adv/(mmadv^7)))/gdpadv;

taxesgdplac=(taulac*(wagddslac*labslac+wagddulac*labddulac)
+tclac*(philac*xs1lac
+xs2lac*philac*p2lac/mmlac
+xs3lac*philac*p3lac/(mmlac^2)
+xs4lac*philac*p4lac/(mmlac^3)
+xs5lac*philac*p5lac/(mmlac^4)
+xs6lac*philac*p6lac/(mmlac^5)
+xs7lac*philac*p7lac/(mmlac^6)
+xs8lac*philac*p8lac/(mmlac^7))
+tclac*((1-philac)*xu1lac
+xu2lac*(1-philac)*p2lac/mmlac
+xu3lac*(1-philac)*p3lac/(mmlac^2)
+xu4lac*(1-philac)*p4lac/(mmlac^3)
+xu5lac*(1-philac)*p5lac/(mmlac^4)
+xu6lac*(1-philac)*p6lac/(mmlac^5)
+xu7lac*(1-philac)*p7lac/(mmlac^6)
+xu8lac*(1-philac)*p8lac/(mmlac^7)))/gdplac;

taxesgdpjap=(taujap*(wagddsjap*labsjap+wagddujap*labddujap)
+tcjap*(phijap*xs1jap
+xs2jap*phijap*p2jap/mmjap
+xs3jap*phijap*p3jap/(mmjap^2)
+xs4jap*phijap*p4jap/(mmjap^3)
+xs5jap*phijap*p5jap/(mmjap^4)
+xs6jap*phijap*p6jap/(mmjap^5)
+xs7jap*phijap*p7jap/(mmjap^6)
+xs8jap*phijap*p8jap/(mmjap^7))
+tcjap*((1-phijap)*xu1jap
+xu2jap*(1-phijap)*p2jap/mmjap
+xu3jap*(1-phijap)*p3jap/(mmjap^2)
+xu4jap*(1-phijap)*p4jap/(mmjap^3)
+xu5jap*(1-phijap)*p5jap/(mmjap^4)
+xu6jap*(1-phijap)*p6jap/(mmjap^5)
+xu7jap*(1-phijap)*p7jap/(mmjap^6)
+xu8jap*(1-phijap)*p8jap/(mmjap^7)))/gdpjap;

taxesgdpeas=(taueas*(wagddseas*labseas+wagddueas*labddueas)
+tceas*(phieas*xs1eas
+xs2eas*phieas*p2eas/mmeas
+xs3eas*phieas*p3eas/(mmeas^2)
+xs4eas*phieas*p4eas/(mmeas^3)
+xs5eas*phieas*p5eas/(mmeas^4)
+xs6eas*phieas*p6eas/(mmeas^5)
+xs7eas*phieas*p7eas/(mmeas^6)
+xs8eas*phieas*p8eas/(mmeas^7))
+tceas*((1-phieas)*xu1eas
+xu2eas*(1-phieas)*p2eas/mmeas
+xu3eas*(1-phieas)*p3eas/(mmeas^2)
+xu4eas*(1-phieas)*p4eas/(mmeas^3)
+xu5eas*(1-phieas)*p5eas/(mmeas^4)
+xu6eas*(1-phieas)*p6eas/(mmeas^5)
+xu7eas*(1-phieas)*p7eas/(mmeas^6)
+xu8eas*(1-phieas)*p8eas/(mmeas^7)))/gdpeas;

taxesgdpnam=(taunam*(wagddsnam*labsnam+wagddunam*labddunam)
+tcnam*(phinam*xs1nam
+xs2nam*phinam*p2nam/mmnam
+xs3nam*phinam*p3nam/(mmnam^2)
+xs4nam*phinam*p4nam/(mmnam^3)
+xs5nam*phinam*p5nam/(mmnam^4)
+xs6nam*phinam*p6nam/(mmnam^5)
+xs7nam*phinam*p7nam/(mmnam^6)
+xs8nam*phinam*p8nam/(mmnam^7))
+tcnam*((1-phinam)*xu1nam
+xu2nam*(1-phinam)*p2nam/mmnam
+xu3nam*(1-phinam)*p3nam/(mmnam^2)
+xu4nam*(1-phinam)*p4nam/(mmnam^3)
+xu5nam*(1-phinam)*p5nam/(mmnam^4)
+xu6nam*(1-phinam)*p6nam/(mmnam^5)
+xu7nam*(1-phinam)*p7nam/(mmnam^6)
+xu8nam*(1-phinam)*p8nam/(mmnam^7)))/gdpnam;

taxesgdpmen=(taumen*(wagddsmen*labsmen+wagddumen*labddumen)
+tcmen*(phimen*xs1men
+xs2men*phimen*p2men/mmmen
+xs3men*phimen*p3men/(mmmen^2)
+xs4men*phimen*p4men/(mmmen^3)
+xs5men*phimen*p5men/(mmmen^4)
+xs6men*phimen*p6men/(mmmen^5)
+xs7men*phimen*p7men/(mmmen^6)
+xs8men*phimen*p8men/(mmmen^7))
+tcmen*((1-phimen)*xu1men
+xu2men*(1-phimen)*p2men/mmmen
+xu3men*(1-phimen)*p3men/(mmmen^2)
+xu4men*(1-phimen)*p4men/(mmmen^3)
+xu5men*(1-phimen)*p5men/(mmmen^4)
+xu6men*(1-phimen)*p6men/(mmmen^5)
+xu7men*(1-phimen)*p7men/(mmmen^6)
+xu8men*(1-phimen)*p8men/(mmmen^7)))/gdpmen;

taxesgdpchi=(tauchi*(wagddschi*labschi+wagdduchi*labdduchi)
+tcchi*(phichi*xs1chi
+xs2chi*phichi*p2chi/mmchi
+xs3chi*phichi*p3chi/(mmchi^2)
+xs4chi*phichi*p4chi/(mmchi^3)
+xs5chi*phichi*p5chi/(mmchi^4)
+xs6chi*phichi*p6chi/(mmchi^5)
+xs7chi*phichi*p7chi/(mmchi^6)
+xs8chi*phichi*p8chi/(mmchi^7))
+tcchi*((1-phichi)*xu1chi
+xu2chi*(1-phichi)*p2chi/mmchi
+xu3chi*(1-phichi)*p3chi/(mmchi^2)
+xu4chi*(1-phichi)*p4chi/(mmchi^3)
+xu5chi*(1-phichi)*p5chi/(mmchi^4)
+xu6chi*(1-phichi)*p6chi/(mmchi^5)
+xu7chi*(1-phichi)*p7chi/(mmchi^6)
+xu8chi*(1-phichi)*p8chi/(mmchi^7)))/gdpchi;

taxesgdpind=(tauind*(wagddsind*labsind+wagdduind*labdduind)
+tcind*(phiind*xs1ind
+xs2ind*phiind*p2ind/mmind
+xs3ind*phiind*p3ind/(mmind^2)
+xs4ind*phiind*p4ind/(mmind^3)
+xs5ind*phiind*p5ind/(mmind^4)
+xs6ind*phiind*p6ind/(mmind^5)
+xs7ind*phiind*p7ind/(mmind^6)
+xs8ind*phiind*p8ind/(mmind^7))
+tcind*((1-phiind)*xu1ind
+xu2ind*(1-phiind)*p2ind/mmind
+xu3ind*(1-phiind)*p3ind/(mmind^2)
+xu4ind*(1-phiind)*p4ind/(mmind^3)
+xu5ind*(1-phiind)*p5ind/(mmind^4)
+xu6ind*(1-phiind)*p6ind/(mmind^5)
+xu7ind*(1-phiind)*p7ind/(mmind^6)
+xu8ind*(1-phiind)*p8ind/(mmind^7)))/gdpind;

taxesgdprus=(taurus*(wagddsrus*labsrus+wagddurus*labddurus)
+tcrus*(phirus*xs1rus
+xs2rus*phirus*p2rus/mmrus
+xs3rus*phirus*p3rus/(mmrus^2)
+xs4rus*phirus*p4rus/(mmrus^3)
+xs5rus*phirus*p5rus/(mmrus^4)
+xs6rus*phirus*p6rus/(mmrus^5)
+xs7rus*phirus*p7rus/(mmrus^6)
+xs8rus*phirus*p8rus/(mmrus^7))
+tcrus*((1-phirus)*xu1rus
+xu2rus*(1-phirus)*p2rus/mmrus
+xu3rus*(1-phirus)*p3rus/(mmrus^2)
+xu4rus*(1-phirus)*p4rus/(mmrus^3)
+xu5rus*(1-phirus)*p5rus/(mmrus^4)
+xu6rus*(1-phirus)*p6rus/(mmrus^5)
+xu7rus*(1-phirus)*p7rus/(mmrus^6)
+xu8rus*(1-phirus)*p8rus/(mmrus^7)))/gdprus;

taxesgdpssa=(taussa*(wagddsssa*labsssa+wagddussa*labddussa)
+tcssa*(phissa*xs1ssa
+xs2ssa*phissa*p2ssa/mmssa
+xs3ssa*phissa*p3ssa/(mmssa^2)
+xs4ssa*phissa*p4ssa/(mmssa^3)
+xs5ssa*phissa*p5ssa/(mmssa^4)
+xs6ssa*phissa*p6ssa/(mmssa^5)
+xs7ssa*phissa*p7ssa/(mmssa^6)
+xs8ssa*phissa*p8ssa/(mmssa^7))
+tcssa*((1-phissa)*xu1ssa
+xu2ssa*(1-phissa)*p2ssa/mmssa
+xu3ssa*(1-phissa)*p3ssa/(mmssa^2)
+xu4ssa*(1-phissa)*p4ssa/(mmssa^3)
+xu5ssa*(1-phissa)*p5ssa/(mmssa^4)
+xu6ssa*(1-phissa)*p6ssa/(mmssa^5)
+xu7ssa*(1-phissa)*p7ssa/(mmssa^6)
+xu8ssa*(1-phissa)*p8ssa/(mmssa^7)))/gdpssa;

// --------- end ----------- growth gdp ------------------------------//

growthadv=(gdpadv/gdpadv-1);
growthnam=(gdpnam/gdpnam-1);
growthjap=(gdpjap/gdpjap-1);
growthssa=(gdpssa/gdpssa-1);
growthlac=(gdplac/gdplac-1);
growthmen=(gdpmen/gdpmen-1);
growthrus=(gdprus/gdprus-1);
growtheas=(gdpeas/gdpeas-1);
growthchi=(gdpchi/gdpchi-1);
growthind=(gdpind/gdpind-1);

// -------- end ---------- ownership ratio ---------------------------//

ownadv=weaadv/kadv;
ownnam=weanam/knam;
ownjap=weajap/kjap;
ownssa=weassa/kssa;
ownlac=wealac/klac;
ownmen=weamen/kmen;
ownrus=wearus/krus;
owneas=weaeas/keas;
ownchi=weachi/kchi;
ownind=weaind/kind;

ownsum=ownadv+ownnam+ownjap+ownssa+ownlac+ownmen+ownrus+owneas+ownchi+ownind;

// -------- end ---------- foreign assets fa ------------------------- //

faadv=weaadv-kadv;
fanam=weanam-knam;
fajap=weajap-kjap;
falac=wealac-klac;
faeas=weaeas-keas;
famen=weamen-kmen;
farus=wearus-krus;
fachi=weachi-kchi;
faind=weaind-kind;
fassa=weassa-kssa;

fasum=faadv+fanam+fajap+falac+faeas+famen+farus+fachi+faind+fassa;

// --------- end --------- courrent account ca ----------------------- //

caadv=faadv-faadv/(ggnam*mmadv);
canam=fanam-fanam/(ggnam*mmnam);
cajap=fajap-fajap/(ggnam*mmjap);
calac=falac-falac/(ggnam*mmlac);
caeas=faeas-faeas/(ggnam*mmeas);
camen=famen-famen/(ggnam*mmmen);
carus=farus-farus/(ggnam*mmrus);
cachi=fachi-fachi/(ggnam*mmchi);
caind=faind-faind/(ggnam*mmind);
cassa=fassa-fassa/(ggnam*mmssa);

casum=caadv+canam+cajap+calac+caeas+camen+carus+cachi+caind+cassa;

// ------------------------- dette intérieure ---------------------------//

zzdebtadv=ddyadv*gdpadv;
zzdebtnam=ddynam*gdpnam;
zzdebtjap=ddyjap*gdpjap;
zzdebtlac=ddylac*gdplac;
zzdebteas=ddyeas*gdpeas;
zzdebtmen=ddymen*gdpmen;
zzdebtrus=ddyrus*gdprus;
zzdebtchi=ddychi*gdpchi;
zzdebtind=ddyind*gdpind;
zzdebtssa=ddyssa*gdpssa;

zzdebtsum=zzdebtadv+zzdebtnam+zzdebtjap+zzdebtlac+zzdebteas+zzdebtmen+zzdebtrus+zzdebtchi+zzdebtind+zzdebtssa;

// ------------------------- zzownership ratio ---------------------------//

zzown2adv=weaadv/(kadv+zzdebtadv);
zzown2nam=weanam/(knam+zzdebtnam);
zzown2jap=weajap/(kjap+zzdebtjap);
zzown2ssa=weassa/(kssa+zzdebtssa);
zzown2lac=wealac/(klac+zzdebtlac);
zzown2men=weamen/(kmen+zzdebtmen);
zzown2rus=wearus/(krus+zzdebtrus);
zzown2eas=weaeas/(keas+zzdebteas);
zzown2chi=weachi/(kchi+zzdebtchi);
zzown2ind=weaind/(kind+zzdebtind);

zzown2sum=zzown2adv+zzown2nam+zzown2jap+zzown2ssa+zzown2lac+zzown2men+zzown2rus+zzown2eas+zzown2chi+zzown2ind;

// ----------------------- foreign assets zzfa2 ------------------------- //

zzfa2adv=weaadv-kadv-zzdebtadv;
zzfa2nam=weanam-knam-zzdebtnam;
zzfa2jap=weajap-kjap-zzdebtjap;
zzfa2lac=wealac-klac-zzdebtlac;
zzfa2eas=weaeas-keas-zzdebteas;
zzfa2men=weamen-kmen-zzdebtmen;
zzfa2rus=wearus-krus-zzdebtrus;
zzfa2chi=weachi-kchi-zzdebtchi;
zzfa2ind=weaind-kind-zzdebtind;
zzfa2ssa=weassa-kssa-zzdebtssa;

zzfa2sum=zzfa2adv+zzfa2nam+zzfa2jap+zzfa2lac+zzfa2eas+zzfa2men+zzfa2rus+zzfa2chi+zzfa2ind+zzfa2ssa;

// ----------------------- courrent account zzca2 ----------------------- //

zzca2adv=zzfa2adv-zzfa2adv/(ggnam*mmadv);
zzca2nam=zzfa2nam-zzfa2nam/(ggnam*mmnam);
zzca2jap=zzfa2jap-zzfa2jap/(ggnam*mmjap);
zzca2lac=zzfa2lac-zzfa2lac/(ggnam*mmlac);
zzca2eas=zzfa2eas-zzfa2eas/(ggnam*mmeas);
zzca2men=zzfa2men-zzfa2men/(ggnam*mmmen);
zzca2rus=zzfa2rus-zzfa2rus/(ggnam*mmrus);
zzca2chi=zzfa2chi-zzfa2chi/(ggnam*mmchi);
zzca2ind=zzfa2ind-zzfa2ind/(ggnam*mmind);
zzca2ssa=zzfa2ssa-zzfa2ssa/(ggnam*mmssa);

zzca2sum=zzca2adv+zzca2nam+zzca2jap+zzca2lac+zzca2eas+zzca2men+zzca2rus+zzca2chi+zzca2ind+zzca2ssa;

// -------------------- new courrent account --------------------------- //

zzvardebtadv=zzdebtadv-zzdebtadv;
zzvardebtnam=zzdebtnam-zzdebtnam;
zzvardebtjap=zzdebtjap-zzdebtjap;
zzvardebtlac=zzdebtlac-zzdebtlac;
zzvardebteas=zzdebteas-zzdebteas;
zzvardebtmen=zzdebtmen-zzdebtmen;
zzvardebtrus=zzdebtrus-zzdebtrus;
zzvardebtind=zzdebtind-zzdebtind;
zzvardebtchi=zzdebtchi-zzdebtchi;
zzvardebtssa=zzdebtssa-zzdebtssa;

zzvarkadv=kadv-kadv;
zzvarknam=knam-knam;
zzvarkjap=kjap-kjap;
zzvarklac=klac-klac;
zzvarkeas=keas-keas;
zzvarkmen=kmen-kmen;
zzvarkrus=krus-krus;
zzvarkind=kind-kind;
zzvarkchi=kchi-kchi;
zzvarkssa=kssa-kssa;

zzvarweaadv=weaadv-weaadv;
zzvarweanam=weanam-weanam;
zzvarweajap=weajap-weajap;
zzvarwealac=wealac-wealac;
zzvarweaeas=weaeas-weaeas;
zzvarweamen=weamen-weamen;
zzvarwearus=wearus-wearus;
zzvarweaind=weaind-weaind;
zzvarweachi=weachi-weachi;
zzvarweassa=weassa-weassa;

zzznewcaadv=zzvarkadv+zzvardebtadv-zzvarweaadv;
zzznewcanam=zzvarknam+zzvardebtnam-zzvarweanam;
zzznewcajap=zzvarkjap+zzvardebtjap-zzvarweajap;
zzznewcalac=zzvarklac+zzvardebtlac-zzvarwealac;
zzznewcaeas=zzvarkeas+zzvardebteas-zzvarweaeas;
zzznewcamen=zzvarkmen+zzvardebtmen-zzvarweamen;
zzznewcarus=zzvarkrus+zzvardebtrus-zzvarwearus;
zzznewcaind=zzvarkind+zzvardebtind-zzvarweaind;
zzznewcachi=zzvarkchi+zzvardebtchi-zzvarweachi;
zzznewcassa=zzvarkssa+zzvardebtssa-zzvarweassa;

zzznewcasum=zzznewcaadv+zzznewcanam+zzznewcajap+zzznewcalac+zzznewcaeas+zzznewcamen+zzznewcarus+zzznewcachi+zzznewcaind+zzznewcassa;

// -------------------------- conspubgdp -------------------------------//

conspubgdpadv=conspubadv*gdpadv;
conspubgdpnam=conspubnam*gdpnam;
conspubgdpjap=conspubjap*gdpjap;
conspubgdpssa=conspubssa*gdpssa;
conspubgdplac=conspublac*gdplac;
conspubgdpmen=conspubmen*gdpmen;
conspubgdprus=conspubrus*gdprus;
conspubgdpeas=conspubeas*gdpeas;
conspubgdpchi=conspubchi*gdpchi;
conspubgdpind=conspubind*gdpind;

end;

steady(solve_algo=3);
//check;


shocks;

// -----------  adv  ------- shocks on population growth ------------- //

var p2adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.737909705     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.583190133     0.700937377     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.503607555     0.603885033     0.724129591
0.931150656     0.931150656     0.931150656     0.931150656     0.931150656
 0.931150656     0.913727623     0.94237471      0.965192065 0.999853812     1
      1       1       1       1       1       1       1     1       1       1
     1       1
;
var p5adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.378770303     0.4541904       0.544628018
0.653073421     0.83978027      0.83978027      0.83978027      0.83978027
0.83978027      0.851911004     0.84525192      0.883235148 0.91763412
0.959286602     0.996323771     1       1       1       1     1       1
1       1       1       1       1
;
var p6adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.246313345     0.295358838     0.354170186
0.42469195      0.509255887     0.654846811     0.654846811     0.654846811
 0.654846811     0.666005536     0.67788039      0.724291031 0.768026054
0.814798973     0.86245523      0.90668286      0.95239868     0.976430451
0.976430451     0.976430451     0.976430451 0.976430451     0.976430451
0.976430451     0.976430451     0.976430451
;
var p7adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.106602647     0.127829184     0.153282313
0.183803626     0.220402291     0.264288419     0.33984571      0.33984571
0.33984571      0.358716145     0.390600311     0.433508694 0.498405133
0.551132592     0.602142653     0.656421314     0.708330598     0.75710297
 0.75710297      0.75710297      0.75710297 0.75710297      0.75710297
0.75710297      0.75710297      0.75710297
;
var p8adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.018188066     0.021809643     0.026152342
0.031359751     0.037604051     0.045091706     0.05407029      0.069528419
 0.069528419     0.076933959     0.093751186     0.117228621 0.151445195
0.207377725     0.240478759     0.283168357     0.330385206     0.379784022
 0.379784022     0.379784022     0.379784022 0.379784022     0.379784022
0.379784022     0.379784022     0.379784022
;

// -----------  nam  ------- shocks on population growth ------------- //

var p2nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.636432622     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.461119686     0.604544715     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.390384311     0.510622679     0.667894464
0.99577746      0.99577746      0.99577746      0.99577746      0.99577746
0.99577746      0.993548831     1       1       1       1       1  1       1
    1       1       1       1       1       1       1       1       1
;
var p5nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.283128152     0.370331622     0.484393762
0.633587041     0.944627824     0.944627824     0.944627824     0.944627824
 0.944627824     0.899145345     0.909478155     0.940986163 0.982585015
0.997363103     1       1       1       1       1       1     1       1
1       1       1       1
;
var p6nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.170478555     0.22298595      0.291665622
0.381498634     0.499000213     0.743969582     0.743969582     0.743969582
 0.743969582     0.740136849     0.729172153     0.76296682 0.803547067
0.855699287     0.878407019     0.895466522     0.962934307     0.98579086
 0.98579086      0.98579086      0.98579086 0.98579086      0.98579086
0.98579086      0.98579086      0.98579086
;
var p7nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.07064341      0.09240158      0.120861267
0.158086538     0.206777191     0.270464566     0.403241131     0.403241131
 0.403241131     0.433041351     0.462067849     0.46980248 0.520100166
0.559656987     0.62013115      0.651387451     0.676592363     0.742740395
 0.742740395     0.742740395     0.742740395 0.742740395     0.742740395
0.742740395     0.742740395     0.742740395
;
var p8nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.014933077     0.019532465     0.025548464
0.033417391     0.043709947     0.057172611     0.074781775     0.111493672
 0.111493672     0.113423494     0.147333468     0.155385172 0.176866312
0.219578588     0.246267784     0.292471423     0.320666371     0.344004851
 0.344004851     0.344004851     0.344004851 0.344004851     0.344004851
0.344004851     0.344004851     0.344004851
;

// -----------  eas  ------- shocks on population growth ------------- //

var p2eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.847752891     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.747812449     0.805276964     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.685932238     0.736930917     0.79172132
0.886155508     0.886155508     0.886155508     0.886155508     0.886155508
 0.886155508     0.898787349     0.893099482     0.899570219 0.882894021
0.888530796     0.901098334     0.946485149     0.956279995     0.959572785
 0.959572785     0.959572785     0.959572785 0.959572785     0.959572785
0.959572785     0.959572785     0.959572785
;
var p5eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.562651715     0.604484556     0.649427645
0.697712229     0.780933289     0.780933289     0.780933289     0.780933289
 0.780933289     0.798962505     0.809848059     0.797124784 0.801036531
0.806735572     0.820748926     0.842726757     0.892376238     0.90729252
 0.90729252      0.90729252      0.90729252 0.90729252      0.90729252
0.90729252      0.90729252      0.90729252
;
var p6eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.385404591     0.414059208     0.444844279
0.477918203     0.513451155     0.574694097     0.574694097     0.574694097
 0.574694097     0.605645259     0.620117767     0.635896669 0.624209316
0.647958539     0.670243806     0.691119296     0.724797607     0.778366337
 0.778366337     0.778366337     0.778366337 0.778366337     0.778366337
0.778366337     0.778366337     0.778366337
;
var p7eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.17026944      0.182928878     0.196529538
0.211141399     0.226839645     0.243705046     0.272773466     0.272773466
 0.272773466     0.304865436     0.325711343     0.329875475 0.366075656
0.371880391     0.402540766     0.439048514     0.464804584     0.508173052
 0.508173052     0.508173052     0.508173052 0.508173052     0.508173052
0.508173052     0.508173052     0.508173052
;
var p8eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.027344764     0.029377832     0.031562057
0.033908678     0.03642977      0.039138303     0.042048214     0.047063601
 0.047063601     0.056094543     0.066319608     0.06907781 0.07933714
0.109037259     0.116158712     0.136518771     0.167499135     0.187129118
 0.187129118     0.187129118     0.187129118 0.187129118     0.187129118
0.187129118     0.187129118     0.187129118
;

// -----------  men  ------- shocks on population growth ------------- //

var p2men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.883787741     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.745430149     0.839506319     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.568707292     0.638996957     0.717974109
0.77084612      0.77084612      0.77084612      0.77084612      0.77084612
0.77084612      0.773019711     0.845076122     0.822999064 0.880943706
0.937650933     0.931637048     0.95338357      0.959676172     0.965180767
 0.965180767     0.965180767     0.965180767 0.965180767     0.965180767
0.965180767     0.965180767     0.965180767
;
var p5men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.419926284     0.471827285     0.530143017
0.595666311     0.639531509     0.639531509     0.639531509     0.639531509
 0.639531509     0.670954811     0.688488346     0.74447423 0.746789147
0.793465719     0.863051592     0.869329901     0.898480595     0.909809812
 0.909809812     0.909809812     0.909809812 0.909809812     0.909809812
0.909809812     0.909809812     0.909809812
;
var p6men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.248344451     0.279038709     0.31352664
0.352277123     0.395816992     0.424965175     0.424965175     0.424965175
 0.424965175     0.480571715     0.486989614     0.520822279 0.58848819
0.591647529     0.652760066     0.724939627     0.744484596     0.782296592
 0.782296592     0.782296592     0.782296592 0.782296592     0.782296592
0.782296592     0.782296592     0.782296592
;
var p7men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.084274197     0.094690109     0.106393381
0.119543124     0.134318117     0.150919233     0.162033009     0.162033009
 0.162033009     0.200992215     0.22419106      0.235235463 0.274344516
0.31606899      0.338776899     0.402562701     0.466169045     0.501229312
 0.501229312     0.501229312     0.501229312 0.501229312     0.501229312
0.501229312     0.501229312     0.501229312
;
var p8men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.007979153     0.00896534      0.010073416
0.011318445     0.012717354     0.014289161     0.016055238     0.017237554
 0.017237554     0.026604546     0.033697155     0.03730206 0.044335534
0.057994768     0.073440898     0.091690057     0.126890144     0.162832053
 0.162832053     0.162832053     0.162832053 0.162832053     0.162832053
0.162832053     0.162832053     0.162832053
;

// -----------  lac  ------- shocks on population growth ------------- //

var p2lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.750142354     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.556601439     0.712557119     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.395251537     0.504826221     0.644778044
0.815592896     0.815592896     0.815592896     0.815592896     0.815592896
 0.815592896     0.821080631     0.831050378     0.839664438 0.853066858
0.859878083     0.87368216      0.892479302     0.909529124     0.918592142
 0.918592142     0.918592142     0.918592142 0.918592142     0.918592142
0.918592142     0.918592142     0.918592142
;
var p5lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.262137139     0.334808821     0.427627108
0.546177198     0.690870674     0.690870674     0.690870674     0.690870674
 0.690870674     0.702579        0.720750021     0.734479715 0.754699174
0.779553018     0.795940704     0.816396223     0.839234545     0.860181373
 0.860181373     0.860181373     0.860181373 0.860181373     0.860181373
0.860181373     0.860181373     0.860181373
;
var p6lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.147394629     0.188256507     0.240446429
0.307104845     0.392242822     0.496155944     0.496155944     0.496155944
 0.496155944     0.513647055     0.539511031     0.563105549 0.593219785
0.626897365     0.663439722     0.689262954     0.716400807     0.743582371
 0.743582371     0.743582371     0.743582371 0.743582371     0.743582371
0.743582371     0.743582371     0.743582371
;
var p7lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.054555101     0.069679288     0.088996318
0.113668564     0.145180641     0.18542874      0.234552594     0.234552594
 0.234552594     0.254514149     0.277390762     0.30012199 0.350864505
0.391447906     0.433844541     0.477045833     0.509310058     0.539856069
 0.539856069     0.539856069     0.539856069 0.539856069     0.539856069
0.539856069     0.539856069     0.539856069
;
var p8lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.007563576     0.009660409     0.012338542
0.015759128     0.020127996     0.025708034     0.032835014     0.041533678
 0.041533678     0.047338822     0.057355301     0.070839207 0.096413557
0.13682587      0.16407859      0.194756849     0.226119347     0.24944583
 0.24944583      0.24944583      0.24944583 0.24944583      0.24944583
0.24944583      0.24944583      0.24944583
;

// -----------  jap  -------- shocks on population growth ------------- //

var p2jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.6281131       0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.415373856     0.596642035     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.289803583     0.415309071     0.595167331
0.899100241     0.899100241     0.899100241     0.899100241     0.899100241
 0.899100241     0.9062793       0.926384484     0.964227367 0.978821038
0.971988622     0.986054929     1       1       1       1     1       1
1       1       1       1       1
;
var p5jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.180504043     0.258675085     0.370699728
0.53123898      0.802525727     0.802525727     0.802525727     0.802525727
 0.802525727     0.813935938     0.847308204     0.878628934 0.924792969
0.941556149     0.942926045     0.961730892     0.979751553     0.99240326
 0.99240326      0.99240326      0.99240326 0.99240326      0.99240326
0.99240326      0.99240326      0.99240326
;
var p6jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.093402302     0.133852118     0.19181957
0.274891038     0.393938338     0.595110041     0.595110041     0.595110041
 0.595110041     0.622826852     0.69076607      0.750909893 0.79019273
0.841586389     0.866017851     0.878679199     0.904566745     0.926086957
 0.926086957     0.926086957     0.926086957 0.926086957     0.926086957
0.926086957     0.926086957     0.926086957
;
var p7jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.031027247     0.044464243     0.063720412
0.091315866     0.130862107     0.187534673     0.283302629     0.283302629
 0.283302629     0.307235477     0.374786664     0.477258366 0.565125593
0.614845084     0.67517323      0.707982452     0.737818452     0.773898233
 0.773898233     0.773898233     0.773898233 0.773898233     0.773898233
0.773898233     0.773898233     0.773898233
;
var p8jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.004218908     0.006045995     0.008664339
0.012416611     0.017793884     0.025499897     0.036543159     0.055204581
 0.055204581     0.056233924     0.080934966     0.130866342 0.207140279
0.289110447     0.339046109     0.401216833     0.438555796     0.483304477
 0.483304477     0.483304477     0.483304477 0.483304477     0.483304477
0.483304477     0.483304477     0.483304477
;

// -----------  ssa  ------- shocks on population growth ------------- //

var p2ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.902353166     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.738544114     0.85714154      0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.498549334     0.57726765      0.668415173
0.702872633     0.702872633     0.702872633     0.702872633     0.702872633
 0.702872633     0.704654295     0.727023117     0.752795833 0.756215419
0.687287383     0.628534924     0.657570715     0.710203709     0.759653392
 0.759653392     0.759653392     0.759653392 0.759653392     0.759653392
0.759653392     0.759653392     0.759653392
;
var p5ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.338863948     0.392368782     0.454321748
0.526056761     0.553175505     0.553175505     0.553175505     0.553175505
 0.553175505     0.565793179     0.579748533     0.609518467 0.633745938
0.626263371     0.576460141     0.539592157     0.576860788     0.634364253
 0.634364253     0.634364253     0.634364253 0.634364253     0.634364253
0.634364253     0.634364253     0.634364253
;
var p6ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.18145258      0.210102987     0.243277143
0.281689324     0.326166585     0.342980794     0.342980794     0.342980794
 0.342980794     0.35753283      0.380463572     0.402405576 0.431578637
0.454884343     0.462154187     0.437680169     0.422007827     0.463109632
 0.463109632     0.463109632     0.463109632 0.463109632     0.463109632
0.463109632     0.463109632     0.463109632
;
var p7ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.05461183      0.063234751     0.073219185
0.084780109     0.098166442     0.113666406     0.119526021     0.119526021
 0.119526021     0.128859076     0.145041609     0.1599109 0.180536364
0.202503617     0.225172051     0.241718923     0.239633855     0.241577047
 0.241577047     0.241577047     0.241577047 0.241577047     0.241577047
0.241577047     0.241577047     0.241577047
;
var p8ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.00494107      0.005721239     0.006624593
0.007670581     0.008881726     0.010284103     0.011907909     0.012521774
 0.012521774     0.01408529      0.017172335     0.02039561 0.025464097
0.03334933      0.040447581     0.049058497     0.057174474     0.059963112
 0.059963112     0.059963112     0.059963112 0.059963112     0.059963112
0.059963112     0.059963112     0.059963112
;

// -----------  rus  ------- shocks on population growth ------------- //

var p2rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.779655487     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.651536269     0.740591522     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.587062227     0.66575932      0.755005947
0.918870933     0.918870933     0.918870933     0.918870933     0.918870933
 0.918870933     0.929225222     0.903051964     0.873271822 0.88541803
0.850482041     0.823926712     0.830798029     0.861320247     0.870982898
 0.870982898     0.870982898     0.870982898 0.870982898     0.870982898
0.870982898     0.870982898     0.870982898
;
var p5rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.457712116     0.519069518     0.588652027
0.667562237     0.812448614     0.812448614     0.812448614     0.812448614
 0.812448614     0.816684995     0.87277968      0.798871918 0.73440638
0.756827525     0.740679123     0.732193732     0.750091443     0.784798781
 0.784798781     0.784798781     0.784798781 0.784798781     0.784798781
0.784798781     0.784798781     0.784798781
;
var p6rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.304849621     0.345715441     0.392059422
0.444615924     0.50421775      0.613652166     0.613652166     0.613652166
 0.613652166     0.633713029     0.632490448     0.677077928 0.573674364
0.531332771     0.578932209     0.579750021     0.588564692     0.615551778
 0.615551778     0.615551778     0.615551778 0.615551778     0.615551778
0.615551778     0.615551778     0.615551778
;
var p7rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.139703701     0.158431316     0.179669412
0.203754526     0.231068306     0.262043565     0.318916978     0.318916978
 0.318916978     0.328884131     0.351561699     0.342250284 0.359622113
0.298718343     0.288762115     0.346311958     0.356432898     0.376706946
 0.376706946     0.376706946     0.376706946 0.376706946     0.376706946
0.376706946     0.376706946     0.376706946
;
var p8rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.027761084     0.031482524     0.035702832
0.040488882     0.045916513     0.052071731     0.059052071     0.071868615
 0.071868615     0.076156935     0.088095187     0.088901255 0.084454347
0.0961499       0.08214894      0.084611578     0.120368673     0.128615306
 0.128615306     0.128615306     0.128615306 0.128615306     0.128615306
0.128615306     0.128615306     0.128615306
;

// -----------  chi  ------- shocks on population growth ------------- //

var p2chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.836132961     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.683547181     0.794239241     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.491653627     0.569947783     0.660710016
0.749794805     0.749794805     0.749794805     0.749794805     0.749794805
 0.749794805     0.767346148     0.837291833     0.897669459 0.924235766
0.927326659     0.928582001     0.930409222     0.934926285     0.943012567
 0.943012567     0.943012567     0.943012567 0.943012567     0.943012567
0.943012567     0.943012567     0.943012567
;
var p5chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.308911372     0.35810445      0.415131357
0.481239604     0.546126056     0.546126056     0.546126056     0.546126056
 0.546126056     0.617668188     0.671554824     0.752348239 0.82106626
0.858361466     0.868551545     0.873828738     0.881553809     0.891094232
 0.891094232     0.891094232     0.891094232 0.891094232     0.891094232
0.891094232     0.891094232     0.891094232
;
var p6chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.16753083      0.19420954      0.225136745
0.260989        0.302550605     0.343344079     0.343344079     0.343344079
0.343344079     0.355765139     0.451398943     0.50847355      0.58936333
0.665696877     0.718604759     0.740631634     0.753983577 0.772421699
0.772421699     0.772421699     0.772421699     0.772421699     0.772421699
 0.772421699     0.772421699     0.772421699
;
var p7chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.054191479     0.062821286     0.072825361
0.08442255      0.097866552     0.113451466     0.128748343     0.128748343
 0.128748343     0.147531918     0.162627249     0.218362137 0.26796558
0.339587166     0.406549886     0.467778875     0.500926956     0.521311169
 0.521311169     0.521311169     0.521311169 0.521311169     0.521311169
0.521311169     0.521311169     0.521311169
;
var p8chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.005700256     0.006608003     0.007660305
0.008880182     0.01029432      0.011933656     0.013834049     0.015699321
 0.015699321     0.018607093     0.02296815      0.033239529 0.045978447
0.072245484     0.103892298     0.137159177     0.176628117     0.202793499
 0.202793499     0.202793499     0.202793499 0.202793499     0.202793499
0.202793499     0.202793499     0.202793499
;

// -----------  ind  -------- shocks on population growth ------------- //

var p2ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.959325838     0.949864173     0.950909591
0.951317554     0.95117407      0.952222059     0.952616763     0.950409741
 0.950032412     0.95982936      0.972921729     0.978023109 0.972715466
0.973405737     0.976937955     0.979571631     0.982535753     0.985385192
 0.984810427     0.984994912     0.985187772 0.985389345     0.985599978
0.985820028     0.986049861     0.986289852
;
var p3ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.877249138     0.911259646     0.903504933
0.905772372     0.907244185     0.906864312     0.908402897     0.908387086
 0.886189877     0.902719134     0.922981791     0.944931726 0.946159968
0.936900814     0.937557748     0.944264423     0.950474173     0.957863989
 0.955282509     0.952976695     0.953299945 0.953637869     0.95399106
0.954360125     0.954745691     0.955148399
;
var p4ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.714826409     0.740820097     0.767759009
0.759391375     0.759391375     0.759391375     0.759391375     0.759391375
 0.759391375     0.764418728     0.788893603     0.816787491 0.855256765
0.881993314     0.890468302     0.909074626     0.928166104     0.941906608
 0.941906608     0.941906608     0.941906608 0.941906608     0.941906608
0.941906608     0.941906608     0.941906608
;
var p5ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.555043112     0.575226498     0.596143825
0.617821782     0.61108828      0.61108828      0.61108828      0.61108828
0.61108828      0.628831302     0.650372481     0.672116155 0.71456104
0.763784067     0.800683037     0.82039103      0.848005054     0.873350269
 0.873350269     0.873350269     0.873350269 0.873350269     0.873350269
0.873350269     0.873350269     0.873350269
;
var p6ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.335815557     0.348027032     0.360682561
0.37379829      0.387390955     0.383168868     0.383168868     0.383168868
 0.383168868     0.407277714     0.438134553     0.45109628 0.48986942
0.541557468     0.601173856     0.647937661     0.680201732     0.717620394
 0.717620394     0.717620394     0.717620394 0.717620394     0.717620394
0.717620394     0.717620394     0.717620394
;
var p7ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.110766231     0.114794094     0.118968424
0.123294549     0.127777987     0.132424459     0.130981195     0.130981195
 0.130981195     0.145233361     0.167534599     0.195433492 0.22324047
0.25954839      0.30501012      0.359393119     0.405151639     0.442768491
 0.442768491     0.442768491     0.442768491 0.442768491     0.442768491
0.442768491     0.442768491     0.442768491
;
var p8ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     0.012399205     0.012850085     0.013317361
0.013801629     0.014303506     0.014823633     0.015362675     0.01519524
0.01519524      0.017061322     0.020650202     0.030745102 0.042858973
0.057869269     0.07338896      0.093001314     0.120708978     0.146184372
 0.146184372     0.146184372     0.146184372 0.146184372     0.146184372
0.146184372     0.146184372     0.146184372
;

//------------------- shocks on debtratio ------------------------------//

var ddynam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.042809925     0.042809925     0.042809925
0.042809925 0.042809925     0.042809925     0.042809925     0.042809925
0.040669429     0.035318189     0.039093472     0.047024413     0.060992156
0.060992156     0.060992156     0.060992156     0.060992156     0.060992156
 0.060992156     0.060992156     0.060992156     0.060992156 0.060992156
0.060992156     0.060992156     0.060992156
;
var ddyadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.044976211     0.044976211     0.044976211
0.044976211 0.044976211     0.044976211     0.044976211     0.044976211
0.044976211     0.041602995     0.038229779     0.051722642     0.066182791
0.066182791     0.066182791     0.066182791     0.066182791     0.066182791
 0.066182791     0.066182791     0.066182791     0.066182791 0.066182791
0.066182791     0.066182791     0.066182791
;
var ddyeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.002290146     0.002290146     0.002290146
0.002290146 0.002290146     0.002290146     0.002290146     0.002290146
0.002290146     0.002290146     0.002290146     0.002290146     0.003489057
0.003489057     0.003489057     0.003489057     0.003489057     0.003489057
 0.003489057     0.003489057     0.003489057     0.003489057 0.003489057
0.003489057     0.003489057     0.003489057
;
var ddymen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.001932557     0.001932557     0.001932557
0.001932557 0.001932557     0.001932557     0.001932557     0.001932557
0.001546046     0.001835929     0.002653999     0.004999042     0.004304415
0.004304415     0.004304415     0.004304415     0.004304415     0.004304415
 0.004304415     0.004304415     0.004304415     0.004304415 0.004304415
0.004304415     0.004304415     0.004304415
;
var ddylac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.001486455     0.001486455     0.001486455
0.001486455 0.001486455     0.001486455     0.001486455     0.001486455
0.001486455     0.001486455     0.004012812     0.003674548     0.004104364
0.004104364     0.004104364     0.004104364     0.004104364     0.004104364
 0.004104364     0.004104364     0.004104364     0.004104364 0.004104364
0.004104364     0.004104364     0.004104364
;
var ddyjap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.068588        0.068588        0.068588
0.068588 0.068588        0.06927388      0.0703027       0.0668733
0.0651586     0.068588        0.0720174       0.085735        0.13396642
0.13396642      0.13396642      0.13396642      0.13396642      0.13396642
0.13396642      0.13396642      0.13396642      0.13396642 0.13396642
0.13396642      0.13396642      0.13396642
;
var ddyssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.001456596     0.001456596     0.001456596
0.001456596 0.001456596     0.001310936     0.001456596     0.001456596
0.001456596     0.001383766     0.002435258     0.004958908     0.003219776
0.003219776     0.003219776     0.003219776     0.003219776     0.003219776
 0.003219776     0.003219776     0.003219776     0.003219776 0.003219776
0.003219776     0.003219776     0.003219776
;
var ddyrus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.003005931     0.003005931     0.003005931
0.003005931 0.003005931     0.003005931     0.003005931     0.003005931
0.003005931     0.003005931     0.003005931     0.002479893     0.003005931
0.003005931     0.003005931     0.003005931     0.003005931     0.003005931
 0.003005931     0.003005931     0.003005931     0.003005931 0.003005931
0.003005931     0.003005931     0.003005931
;
var ddychi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.00099159      0.00099159      0.00099159
0.00099159 0.00099159      0.00099159      0.00099159      0.001090749
0.00099159     0.001090749     0.00098226      0.001930948     0.001833757
0.001833757     0.001833757     0.001833757     0.001833757     0.001833757
 0.001833757     0.001833757     0.001833757     0.001833757 0.001833757
0.001833757     0.001833757     0.001833757
;
var ddyind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.000955737     0.000955737     0.000955737
0.000955737 0.000955737     0.000955737     0.000955737     0.000955737
0.000955737     0.000955737     0.001399205     0.002562997     0.002531276
0.002531276     0.002531276     0.002531276     0.002531276     0.002531276
 0.002531276     0.002531276     0.002531276     0.002531276 0.002531276
0.002531276     0.002531276     0.002531276
;

//--------------------- shocks on conspublique -------------------------//

var conspubnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.170280977     0.170280977     0.170280977
0.170280977     0.170280977     0.170280977     0.170280977     0.170280977
 0.164321143 0.169202777     0.15708339      0.159702597     0.149242875
0.149242875     0.149242875     0.149242875     0.149242875     0.149242875
0.149242875     0.149242875     0.149242875     0.149242875     0.149242875
 0.149242875     0.149242875     0.149242875
;
var conspubadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.142865046     0.142865046     0.142865046
0.142865046     0.142865046     0.142865046     0.142865046     0.142865046
 0.142865046 0.163146515     0.197673806     0.198800038     0.195747583
0.195747583     0.195747583     0.195747583     0.195747583     0.195747583
0.171279136     0.171279136     0.171279136     0.171279136     0.171279136
 0.171279136     0.171279136     0.171279136
;
var conspubeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.104169516     0.104169516     0.104169516
0.104169516     0.104169516     0.104169516     0.104169516     0.105732058
 0.104169516 0.106773754     0.117547985     0.136812642     0.164212094
0.164212094     0.164212094     0.164212094     0.164212094     0.164212094
0.164212094     0.164212094     0.164212094     0.164212094     0.164212094
 0.164212094     0.164212094     0.164212094
;
var conspubmen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.167804208     0.167804208     0.167804208
0.167804208     0.167804208     0.167804208     0.167804208     0.167804208
 0.167804208 0.161295015     0.161988101     0.151067262     0.154196071
0.154196071     0.154196071     0.154196071     0.154196071     0.154196071
0.154196071     0.154196071     0.154196071     0.154196071     0.154196071
 0.154196071     0.154196071     0.154196071
;
var conspublac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.101727665     0.101727665     0.101727665
0.101727665     0.101727665     0.101727665     0.101727665     0.101727665
 0.101727665 0.102440887     0.106998074     0.117975453     0.150334205
0.150334205     0.150334205     0.150334205     0.150334205     0.150334205
0.150334205     0.150334205     0.150334205     0.150334205     0.150334205
 0.150334205     0.150334205     0.150334205
;
var conspubjap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.108820692     0.108820692     0.108820692
0.108820692     0.108820692     0.108820692     0.108820692     0.110453003
 0.108820692 0.104568715     0.13509521      0.129785657     0.164831366
0.164831366     0.164831366     0.164831366     0.164831366     0.164831366
0.164831366     0.164831366     0.164831366     0.164831366     0.164831366
 0.164831366     0.164831366     0.164831366
;
var conspubssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.100729755     0.100729755     0.100729755
0.100729755     0.100729755     0.100729755     0.100729755     0.100729755
 0.100729755 0.123849429     0.139029694     0.148638944     0.150356571
0.150356571     0.150356571     0.150356571     0.150356571     0.150356571
0.150356571     0.150356571     0.150356571     0.150356571     0.150356571
 0.150356571     0.150356571     0.150356571
;
var conspubrus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.126599505     0.126599505     0.126599505
0.126599505     0.126599505     0.126599505     0.126599505     0.126599505
 0.126599505 0.126599505     0.124067515     0.125006736     0.162501299
0.162501299     0.162501299     0.162501299     0.162501299     0.162501299
0.162501299     0.162501299     0.162501299     0.162501299     0.162501299
 0.162501299     0.162501299     0.162501299
;
var conspubchi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.077531613     0.077531613     0.077531613
0.077531613     0.077531613     0.077531613     0.077531613     0.077531613
 0.077531613 0.089595279     0.108449084     0.114667849     0.125176012
0.125176012     0.125176012     0.125176012     0.125176012     0.125176012
0.125176012     0.125176012     0.125176012     0.125176012     0.125176012
 0.125176012     0.125176012     0.125176012
;
var conspubind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
22 23 24 25 26:49; values 0.080276755     0.080276755     0.080276755
0.080276755     0.080276755     0.080276755     0.080276755     0.080276755
 0.080276755 0.097866454     0.099278094     0.106501051     0.110200901
0.110200901     0.110200901     0.110200901     0.110200901     0.110200901
0.110200901     0.110200901     0.110200901     0.110200901     0.110200901
 0.110200901     0.110200901     0.110200901
;

//------------------ shocks ---------- consumption tax -----------------//

var tcnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.111111111     0.111111111     0.111111111
0.111111111 0.111111111     0.111111111     0.125888889     0.14
0.173333333 0.191111111     0.184166667     0.204166667     0.2     0.2
0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2
 0.2  0.2
;
var tcadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.111111111     0.111111111     0.111111111
0.111111111 0.111111111     0.113888889     0.122222222     0.12
0.144444444 0.177777778     0.188888889     0.194444444     0.2     0.2
0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2
 0.2  0.2
;
var tceas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.073333333
0.079444444     0.097777778     0.109083333     0.106944444     0.11    0.11
 0.11  0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11
0.11    0.11
;
var tcmen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.073333333
0.085402778     0.099733333     0.119472222     0.106944444     0.11    0.11
 0.11  0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11
0.11    0.11
;
var tclac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.073333333
0.083416667     0.088   0.114277778     0.088229167     0.11    0.11    0.11
 0.11  0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11
0.11
;
var tcjap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.111111111     0.111111111     0.111111111
0.111111111 0.111111111     0.111111111     0.122222222     0.133333333
0.144444444     0.177777778     0.179444444     0.170138889     0.2     0.2
 0.2  0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2 0.2
   0.2
;
var tcssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.073333333
0.079444444     0.129555556     0.155833333     0.139027778     0.11    0.11
 0.11  0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11
0.11    0.11
;
var tcrus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.0715
0.079444444 0.097777778     0.103888889     0.157743056     0.11    0.11
0.11 0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11 0.11
   0.11
;
var tcchi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.073333333
0.079444444     0.095333333     0.064930556     0.052402778     0.11    0.11
 0.11  0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11
0.11    0.11
;
var tcind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49; values      0.061111111     0.061111111     0.061111111
0.061111111 0.061111111     0.061111111     0.067222222     0.073333333
0.079444444     0.1012  0.101291667     0.106944444     0.11    0.11    0.11
 0.11  0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11
0.11
;

//--------------- shocks replacement rate rep --------------------------//

var repnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.2075  0.2075  0.2075  0.2075  0.2075  0.2075
0.2075  0.2075 0.259375        0.31125 0.363125        0.415   0.415   0.404
0.393   0.382   0.371   0.36 0.36    0.36    0.36    0.36    0.36    0.36
0.36    0.36
;
var repadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.2125  0.2125  0.2125  0.2125  0.2125  0.2125
0.2125  0.2125 0.265625        0.31875 0.371875        0.425   0.425   0.414
0.403   0.392   0.381   0.37 0.37    0.37    0.37    0.37    0.37    0.37
0.37    0.37
;
var repeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.2125  0.2125  0.2125  0.2125  0.2125  0.2125
0.2125  0.2125 0.265625        0.31875 0.371875        0.425   0.425   0.425
0.425   0.425   0.425   0.425 0.425   0.425   0.425   0.425   0.425   0.425
0.425   0.425
;
var repmen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.2     0.2     0.2     0.2     0.2     0.2     0.2
    0.2 0.25    0.3     0.35    0.4     0.4     0.4     0.4     0.4     0.4
 0.4     0.4     0.4     0.4     0.4     0.4     0.4     0.4     0.4
;
var replac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.125   0.125   0.125   0.125   0.125   0.125
0.125   0.125 0.15625 0.1875  0.21875 0.25    0.25    0.25    0.25    0.25
0.25    0.25    0.25    0.25 0.25    0.25    0.25    0.25    0.25    0.25
;
var repjap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.1375  0.1375  0.1375  0.1375  0.1375  0.1375
0.1375  0.1375 0.171875        0.20625 0.240625        0.275   0.275   0.264
0.253   0.242   0.231   0.22 0.22    0.22    0.22    0.22    0.22    0.22
0.22    0.22
;
var repssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.0375  0.0375  0.0375  0.0375  0.0375  0.0375
0.0375  0.0375 0.046875        0.05625 0.065625        0.075   0.075   0.075
0.075   0.075   0.075   0.075 0.075   0.075   0.075   0.075   0.075   0.075
0.075   0.075
;
var reprus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.225   0.225   0.225   0.225   0.225   0.225
0.225   0.225 0.28125 0.3375  0.39375 0.45    0.45    0.45    0.45    0.45
0.45    0.45    0.45    0.45 0.45    0.45    0.45    0.45    0.45    0.45
;
var repchi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.15    0.15    0.15    0.15    0.15    0.15
0.15    0.15 0.1875  0.225   0.2625  0.3     0.3     0.3     0.3     0.3
0.3     0.3     0.3     0.3     0.3     0.3     0.3     0.3     0.3     0.3

;
var repind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.375   0.375   0.375   0.375   0.375   0.375
0.375   0.375 0.46875 0.5625  0.65625 0.75    0.75    0.75    0.75    0.75
0.75    0.75    0.75    0.75 0.75    0.75    0.75    0.75    0.75    0.75
;

//-------------------------- shocks on mm  -------------------------//

var mmnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.304725656     1.48659229      1.244390949
1.138966681     1.256260388     1.160434015     0.952488625     1.090173354
 1.431325881     1.177891962     0.917263535     1.047459359 1.112045976
0.987444341     1.034684953     1.01409443      1.019806974     1.04964396
 1.04962793      1.049611265     1.049593949 1.049575966     1.049557301
1.049537937     1.049517862     1
;
var mmadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.196116518     1.311395468     1.14941639
1.219559241     1.09070021      0.936104331     1.097953935     1.012073365
 1.113627525     1.140827264     1.018353653     0.873889021 0.955996284
0.941345278     0.961561539     1.001355938     1.003215022     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.071660001     1.182556046     1.209118711
1.379305901     1.163357608     0.937726165     1.320726152     0.900969122
 1.236137585     0.946576569     1.006700167     1.087146758 0.845747172
0.735526865     0.935010124     0.917914682     0.88994537     0.999660915
0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmmen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.120782786     1.21610902      1.188976026
1.255549969     1.215078874     1.19943558      1.292437927     1.170745012
 1.476581934     1.455614356     1.315081885     1.359042294 1.109713323
0.998428051     1.10684688      1.001092799     0.970196439     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmlac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.274030412     1.425960203     1.384041576
1.295090486     1.240877839     1.242174864     1.254134084     1.227105175
 1.363684893     1.361237476     1.230101103     1.173452937 1.065258834
1.019818583     0.988949924     0.966788036     0.945854919     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmjap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.42948343      1.561702179     1.222408219
1.307451543     1.217945727     1.185235075     1.38495464      1.096795003
 1.100373222     0.805927297     1.186612228     0.867508904 0.784828611
0.948653874     0.932539134     0.909851286     0.996252361     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.154996155     1.19012962      1.196092464
1.198904455     1.221043173     1.264114244     1.223842194     1.239098236
 1.287279154     1.295817985     1.236354924     1.384675549 1.354822468
1.265788584     1.246078701     1.137262553     1.088876736     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmrus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.131213491     1.24170553      1.146590632
1.319262583     1.356505614     0.883044864     1.55957394      0.900309783
 1.14168977      1.169268662     0.841395805     1.133901052 0.978563665
0.713515635     1.05313342      0.887962562     0.889089681     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmchi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.156344609     1.526697904     1.223293886
1.029885905     1.101027645     1.136643323     1.187590173     1.109346786
 1.466024095     1.268604948     1.278139903     0.85168693 1.089554455
0.840292383     0.973047699     0.974323038     0.897926599     0.999660915
 0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;
var mmind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26:49;  values     1.033769286     1.164773996     1.191858705
1.177044004     1.266885477     1.218714623     1.249482588     1.169656384
 1.25292936      1.377121195     1.259398912     1.199135153 1.175703687
1.047073439     1.023387394     0.982640009     0.95042704     0.999660915
0.999645647     0.999629776     0.999613285 0.999596158     0.999578382
0.99955994      0.999540821     1
;

//------------------------- shocks on nn  --------------------------//

var nnadvnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   2.737660895     2.509770081     2.213990435
2.045013984 2.189717877     1.901139104     1.533619772     1.767836191
1.641188455     1.276908817     1.236728359     1.373026173     1.14550745
0.984762221     0.938788373     0.872442177     0.861483042     0.847466972
 0.807111402     0.768677526     0.732073834     0.697213175 0.664012548
0.632392903     0.602278955     0.573599005
;
var nneasnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   0.776607208     0.637880368     0.507421766
0.493038906 0.597077583     0.552922591     0.446806948     0.619545059
0.512020374     0.442196733     0.355357774     0.390006489     0.404783524
0.307851049     0.229311879     0.207221462     0.187567959     0.163683168
 0.155888732     0.148465459     0.141395675     0.134662548 0.128250045
0.1221429       0.116326572     0.110787211
;
var nnmennam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   0.643436259     0.552723309     0.45215612
0.432020811 0.476241952     0.460630249     0.476111784     0.646039136
0.693786078     0.715722396     0.884474832     1.26807268      1.645280449
1.641829271     1.660091947     1.775871571     1.753103251     1.667820062
 1.588400059     1.512761961     1.440725677     1.372119693 1.30678066
1.244553009     1.18528858      1.128846267
;
var nnlacnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   0.628185755     0.613406928     0.58838854
0.654419902 0.744124479     0.735012888     0.786787118     1.035956247
1.166078098     1.110972077     1.283901136     1.721782392     1.928886871
1.847732758     1.908312324     1.823961315     1.738875518     1.612779677
 1.535980645     1.462838709     1.393179723     1.326837832 1.263655078
1.203481026     1.146172406     1.091592768
;
var nnjapnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   0.316280356     0.346523061     0.364031095
0.357600321 0.410499358     0.397979546     0.406485256     0.591045002
0.594634975     0.457142857     0.312782428     0.404629029     0.335114945
0.236508024     0.227217114     0.204785863     0.18373504      0.179491288
 0.170944084     0.16280389      0.155051323     0.147667927 0.140636121
0.133939163     0.127561107     0.121486769
;
var nnssanam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   1.078307209     0.954561347     0.764198591
0.734537789 0.773192615     0.751517419     0.818662554     1.051890542
1.195585739     1.075263585     1.182914849     1.594419208     2.107722151
2.567869843     3.291709921     3.964230475     4.44571111      4.746811432
 4.520772792     4.305497897     4.100474188     3.905213512 3.719250964
3.542143775     3.373470262     3.212828821
;
var nnrusnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   1.263917419     1.095832238     0.91531549
0.843378174 0.97688307      1.054834955     0.802688113     1.314295447
1.085398981     0.865762948     0.859424732     0.788340904     0.85339882
0.750962725     0.542636809     0.552312041     0.483616122     0.421626949
 0.401549476     0.382428072     0.364217211     0.346873535 0.330355747
0.314624521     0.299642401     0.285373715
;
var nnchinam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   4.749303854     4.209185189     4.322741514
4.249454939 3.842477413     3.367672736     3.298630238     4.112826918
4.18516129     4.286618005     4.616743291     6.433095393     5.230735895
5.124942421     4.361207919     4.101406247     3.940554726     3.469606497
 3.30438714      3.147035371     2.997176544     2.854453851 2.718527478
2.589073788     2.46578456      2.348366248
;
var nnindnam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values   4.054790781     3.212719968     2.517228631
2.410963258 2.491565288     2.512638231     2.638830744     3.461640362
3.714023769     3.251118063     3.801013792     5.218775687     5.974472739
6.316474122     6.697908943     6.624775552     6.41929323      5.98257319
5.697688752     5.42637024      5.167971657     4.921877769 4.687502637
4.464288226     4.251703072     4.049241021
;

//--------------------- shocks on skilled proportions ------------------//

var phinam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.093717733     0.093717733     0.093717733
0.093717733 0.117147166     0.146433958     0.183042448     0.228803059
0.302113752     0.469195452     0.5     0.525   0.55    0.55    0.55    0.55
 0.55  0.55    0.55    0.55    0.55    0.55    0.55    0.55    0.55    0.55
;
var phiadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.016504496     0.016504496     0.016504496
0.016504496 0.02063062      0.025788275     0.032235344     0.04029418
0.108049582     0.168556391     0.24    0.278653805     0.3     0.3     0.3
 0.3  0.3     0.3     0.3     0.3     0.3     0.3     0.3     0.3     0.3 0.3

;
var phieas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.015494457     0.015494457     0.015494457
0.015494457 0.019368072     0.02421009      0.030262612     0.037828265
0.091359618     0.082126471     0.138541498     0.18584032      0.2     0.2
 0.2  0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2     0.2 0.2
   0.2
;
var phimen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.004347393     0.004347393     0.004347393
0.004347393 0.005434242     0.006792802     0.008491003     0.010613753
0.033885697     0.071755788     0.082896999     0.1332004       0.15    0.15
 0.15  0.15    0.15    0.15    0.15    0.15    0.15    0.15    0.15    0.15
0.15    0.15
;
var philac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.009983981     0.009983981     0.009983981
0.009983981 0.012479976     0.01559997      0.019499963     0.024374954
0.03307769     0.095009906     0.139989615     0.15254428      0.15    0.15
0.15  0.15    0.15    0.15    0.15    0.15    0.15    0.15    0.15    0.15
0.15    0.15
;
var phijap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.03315712      0.03315712      0.03315712
0.03315712 0.0414464       0.051807999     0.064759999     0.080949999
0.04 0.27    0.3     0.32    0.35    0.35    0.35    0.35    0.35    0.35 0.35
   0.35    0.35    0.35    0.35    0.35    0.35    0.35
;
var phissa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.002270809     0.002270809     0.002270809
0.002270809 0.002838512     0.00354814      0.004435175     0.005543968
0.014673323     0.001551918     0.036498519     0.039281761     0.05    0.05
 0.05  0.05    0.05    0.05    0.05    0.05    0.05    0.05    0.05    0.05
0.05    0.05
;
var phirus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.02242942      0.02242942      0.02242942
0.02242942 0.028036775     0.035045969     0.043807461     0.054759326
0.143437283     0.119411389     0.206495405     0.326101236     0.25    0.25
 0.25  0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25
;
var phichi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.006617984     0.006617984     0.006617984
0.006617984 0.00827248      0.0103406       0.01292575      0.016157188
0.013309886     0.022024493     0.061736833     0.071152581     0.1     0.1
 0.1  0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1 0.1
   0.1
;
var phiind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values     0.002329718     0.002329718     0.002329718
0.002329718 0.002912148     0.003640185     0.004550231     0.005687789
0.024016439     0.034087044     0.053144364     0.055489056     0.1     0.1
 0.1  0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1 0.1
   0.1
;

// ------------------------- shocks on psi ------------------------- //

var psinam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psiadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psieas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psimen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psilac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psijap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psissa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psirus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psichi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;
var psiind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:49; values     0.25    0.25    0.25    0.25    0.25    0.25
0.25    0.25 0.4     0.55    0.7     0.95    1       1       1       1       1
      1     1       1       1       1       1       1       1       1
;

// ---------------------- shocks on lams & lamu 5 --------------------- //

var lams5nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.7     0.7     0.7     0.7     0.7     0.7     0.7
  0.7 0.7     0.7     0.7     0.7     0.7     0.72    0.74    0.76    0.78
0.8     0.8     0.8     0.8     0.8     0.8     0.8     0.8     0.8
;
var lams5adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.7     0.7     0.7     0.7     0.7     0.7     0.7
  0.7 0.7     0.7     0.7     0.7     0.7     0.72    0.74    0.76    0.78
0.8     0.8     0.8     0.8     0.8     0.8     0.8     0.8     0.8
;
var lams5eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lams5men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lams5lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lams5jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.7     0.7     0.7     0.7     0.7     0.7     0.7
  0.7 0.7     0.7     0.7     0.7     0.7     0.72    0.74    0.76    0.78
0.8     0.8     0.8     0.8     0.8     0.8     0.8     0.8     0.8
;
var lams5ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lams5rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lams5chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lams5ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;

var lamu5nam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.52    0.54    0.56    0.58
0.6     0.6     0.6     0.6     0.6     0.6     0.6     0.6     0.6
;
var lamu5adv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.52    0.54    0.56    0.58
0.6     0.6     0.6     0.6     0.6     0.6     0.6     0.6     0.6
;
var lamu5eas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lamu5men; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lamu5lac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lamu5jap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.52    0.54    0.56    0.58
0.6     0.6     0.6     0.6     0.6     0.6     0.6     0.6     0.6
;
var lamu5ssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lamu5rus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lamu5chi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;
var lamu5ind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26:33; values   0.5     0.5     0.5     0.5     0.5     0.5     0.5
  0.5 0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5     0.5
;

// ------------- shocks on skill premium -------------------- //

var aaadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.051225354     0.05132346
0.051006191     0.051353981     0.052775478     0.057273857     0.079390417
 0.11020285      0.15364 0.22563639      0.29163103      0.40506609
0.48132536      0.51834073     0.54216328      0.55330774      0.55529325
0.55467315 0.5539995       0.55362385      0.55364874      0.55368488
0.55368681     0.55368882      0.55369091      0.55369309      0.55368663
0.55368081      0.55367564      0.55367118      0.55367118      0.55367118
0.55367118
;
var aachi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.034984369     0.034966913
0.033816958     0.034057256     0.03591725     0.039132271     0.054480762
0.077313776     0.092083555 0.10377004      0.13867089      0.21696072
0.28499411      0.32716509     0.35744666      0.37325825      0.38024337
0.38139566 0.38002483      0.37918308      0.37847744      0.37816992
0.37817422     0.3781787       0.37818338      0.37818825      0.3781792
0.3781727       0.37816878      0.3781675       0.3781675       0.3781675
0.3781675
;
var aaeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.059252714     0.05935244
0.058860925     0.058408217     0.059021848     0.064653089     0.090907789
 0.12462666      0.17224715 0.22292505      0.26067335      0.35630235
0.42305344      0.46490446     0.49542013      0.50526813      0.50658514
0.50583184 0.50340406      0.50198882      0.50098796      0.50063023
0.50063475     0.50063947      0.50064439      0.50064952      0.50063995
0.50063305      0.50062886      0.50062749      0.50062749      0.50062749
0.50062749
;
var aaind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.017332253     0.017312142
0.017084583     0.016881181     0.017268646     0.018659045     0.026406117
 0.038010594     0.064958282 0.10868949      0.14676483      0.21287982
0.27260998      0.3183419     0.34757641      0.36768362      0.37712164
0.3790948 0.37893313      0.37858845      0.37824461      0.37810743
0.37811174     0.37811624      0.37812094      0.37812583      0.37811676
0.37811026      0.37810634      0.37810508      0.37810508      0.37810508
0.37810508
;
var aajap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.099351354     0.099571621
0.098899921     0.10072341      0.10436869     0.11382969      0.15516842
0.21066256      0.23623332 0.30789837      0.39381032      0.49854734
0.5826052       0.63034986     0.64146624      0.65141946      0.65562445
0.65322413 0.65250178      0.65168691      0.65081809      0.65078928
0.65079106     0.65079292      0.65079485      0.65079686      0.65079092
0.65078557      0.65078083      0.65077675      0.65077675      0.65077675
0.65077675
;
var aalac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.044266378     0.044264741
0.043659876     0.043463734     0.04521925     0.04980653      0.069691246
0.098491037     0.12397274 0.182333        0.25261875      0.35438321
0.41285166      0.43601899     0.44592465      0.44920526      0.45073869
0.45212234 0.45155233      0.45096024      0.45052278      0.45036394
0.45036854     0.45037334      0.45037834      0.45038356      0.45037396
0.45036709      0.450363        0.45036167      0.45036167      0.45036167
0.45036167
;
var aamen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.024427944     0.02439484
0.02414939      0.024047192     0.024566425     0.026780679     0.037767165
 0.053852722     0.083835446 0.14601966      0.19664869      0.28918809
0.35992176      0.39748083     0.42214749      0.43253881      0.43678283
0.43853655 0.43894939      0.43904878      0.43887068      0.43878618
0.43879062     0.43879526      0.43880009      0.43880513      0.43879572
0.43878893      0.43878479      0.43878344      0.43878344      0.43878344
0.43878344
;
var aanam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.19824931      0.19895977
0.1974446       0.19924844      0.2071317     0.22244286      0.29182085
0.38102177      0.44690783 0.52990884      0.59318295      0.68516085
0.73293504      0.75266234     0.76233989      0.76663314      0.76884628
0.76943714 0.76889088      0.76868624      0.76840454      0.76821317
0.76821483     0.76821657      0.76821837      0.76822025      0.7687914
0.76926394      0.76966027      0.7699911       0.7699911       0.7699911
0.7699911
;
var aarus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.080930759     0.081167691
0.080583015     0.080918        0.082213217     0.08850542      0.12481007
 0.16742038      0.23509727 0.30135319      0.34417475      0.46859933
0.53633178      0.56438173     0.58681749      0.5810137       0.575027
0.5752283 0.57320177      0.57190822      0.5708485       0.57051859
0.57052331     0.57052824      0.57053338      0.57053874      0.57052905
0.57052221      0.57051824      0.57051699      0.57051699      0.57051699
0.57051699
;
var aassa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
24 25 26 27 28 29 30 31 32 33:49; values 0.01793145      0.017867044
0.017750796     0.017663189     0.018139919     0.01974299      0.027845297
 0.040208054     0.059667278 0.069935812     0.091112601     0.15861224
0.20883917      0.24041485     0.25699231      0.26298821      0.26758372
0.27048431 0.27339569      0.2748902       0.27552723      0.27569453
0.27569871     0.27570308      0.27570763      0.27571238      0.27570415
0.2756985       0.27569545      0.27569454      0.27569454      0.27569454
0.27569454
;

// ---------------------- shocks on unemployment --------------------- //

var etaadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.028799361     0.036368878
    0.046041377     0.058505114   0.074709126     0.096088731     0.028706397
   0.018773315 0.020291803     0.024453205     0.036996943     0.042719526
0.035676528     0.034978042     0.033746848     0.032050714     0.030060795
0.028012313     0.027993451     0.027982935     0.027983632     0.027984644
 0.027984698     0.027984754     0.027984812     0.027984873 0.027984693
0.02798453      0.027984385     0.02798426      0.02798426     0.02798426
0.02798426
;
var etachi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.003354953     0.004203341
    0.005267076     0.006608319   0.008303857     0.010452273     0.007555024
   0.006241689 0.00690921      0.008083529     0.011334506     0.005049266
0.006097755     0.005999323     0.005833913     0.005582555     0.005279234
0.004944844     0.004938227     0.004934173     0.00493078      0.004929302
 0.004929323     0.004929345     0.004929367     0.00492939 0.004929347
0.004929316     0.004929297     0.004929291     0.004929291     0.004929291
 0.004929291
;
var etaeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.00520133      0.006521829
    0.008181477     0.010273663   0.012923514     0.016323428     0.010945181
   0.008797024 0.00989863      0.011915085     0.017043906     0.012317887
0.034692213     0.034043042     0.032965845     0.03118168      0.029144204
0.027078489     0.027011989     0.026973301     0.026945975     0.026936215
 0.026936339     0.026936467     0.026936602     0.026936741 0.02693648
0.026936292     0.026936178     0.026936141     0.026936141     0.026936141
 0.026936141
;
var etaind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.006033912     0.007573918
    0.009516546     0.011973574   0.015092986     0.019071367     0.013748483
   0.01133147 0.012662659     0.015156604     0.021548579     0.012458576
0.014517858     0.014269006     0.013812754     0.013229158     0.012503447
0.011688336     0.011686516     0.011682638     0.011678772     0.01167723
0.011677279     0.01167733      0.011677382     0.011677437 0.011677335
0.011677262     0.011677218     0.011677204     0.011677204     0.011677204
 0.011677204
;
var etajap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.001886256     0.002360675
    0.002953516     0.003702426   0.004648767     0.005863276     0.003983901
   0.003257032 0.003601951     0.004412352     0.006593738     0.005773679
0.01327615     0.013200758     0.012611396     0.011996621     0.011308294
0.010549669     0.010542213     0.010533802     0.010524833     0.010524536
 0.010524554     0.010524573     0.010524593     0.010524614 0.010524552
0.010524497     0.010524448     0.010524406     0.010524406     0.010524406
 0.010524406
;
var etalac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.004394963     0.005509337
    0.006909461     0.00867478   0.010912652     0.013765008     0.009695453
  0.007941499 0.008835512     0.010666661     0.015565321     0.014139132
0.021631624     0.020901287     0.019880911     0.018728509     0.017551646
0.016382597     0.016373121     0.016363286     0.016356025     0.01635339
0.016353466     0.016353546     0.016353629     0.016353716 0.016353556
0.016353442     0.016353374     0.016353353     0.016353353     0.016353353
 0.016353353
;
var etamen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.007953811     0.009989573
    0.01256076      0.015818749   0.0199646       0.025273636     0.016789244
   0.013371911 0.014868411     0.017974314     0.026040068     0.026972246
0.033267444     0.032443542     0.031198558     0.029506931     0.027646662
0.025741205     0.025751785     0.025754333     0.025749767     0.025747602
 0.025747715     0.025747834     0.025747958     0.025748087 0.025747846
0.025747672     0.025747566     0.025747531     0.025747531     0.025747531
 0.025747531
;
var etanam; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.007346759     0.009223341
    0.011567458     0.014567017   0.018470765     0.023624088     0.016533262
   0.014108394 0.016430622     0.020746723     0.030731861     0.020226198
0.013604296     0.013115643     0.012499044     0.011821164     0.011122317
0.010410149     0.010405504     0.010403763     0.010401365     0.010399736
 0.01039975      0.010399765     0.01039978      0.010399796 0.010404657
0.010408676     0.010412045     0.010414856     0.010414856     0.010414856
 0.010414856
;
var etarus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.004312659     0.005407054
    0.006780961     0.008516164   0.010713515     0.013536352     0.009919526
   0.008325402 0.009629154     0.011860082     0.017050648     0.018571051
0.022818989     0.022183581     0.021383213     0.019969589     0.01858101
0.017330581     0.017294668     0.01727177      0.017253026     0.017247193
 0.017247277     0.017247364     0.017247455     0.017247549 0.017247378
0.017247257     0.017247187     0.017247165     0.017247165     0.017247165
 0.017247165
;
var etassa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
23 24 25 26 27 28 29 30 31 32 33:49; values        0.011768195     0.014846073
    0.018773662     0.023813134   0.030328319     0.038833064     0.029739879
   0.025465414 0.029129089     0.035364787     0.051469669     0.047591205
0.051033983     0.048579098     0.045695962     0.042544434     0.03944654
0.03639647      0.036478447     0.036520781     0.036538878     0.036543636
 0.036543755     0.036543879     0.036544008     0.036544143 0.036543909
0.036543748     0.036543662     0.036543636     0.036543636     0.036543636
 0.036543636
;

// ------------------------- shocks on ratiogdp ----------------------- //

var tfpadv; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.54867406
 0.55016499      0.55121556      0.55257325 0.54694673      0.53305848
0.49862037      0.47456139      0.49338723     0.54622994      0.62600955
0.78216742      0.86190594
;
var tfpchi; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.038072378
 0.03827747      0.037949004     0.037950494 0.03742763      0.036270315
0.033675815     0.031491077     0.029989213     0.028510546     0.03180544
 0.058644165     0.12008983
;
var tfpeas; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.29634396
 0.2973081       0.29611524      0.29469976 0.28529983      0.27385948
0.25642445      0.24087838      0.25212456     0.26718522      0.2965297
0.30392631      0.37457905
;
var tfpind; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.044222969
 0.04454361      0.044561181     0.044528604 0.043284816     0.041044388
0.036578589     0.032644145     0.032236979     0.033713483     0.037997868
 0.052396567     0.074057072
;
var tfpjap; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.57740571
 0.57920045      0.58100637      0.58220513 0.57938844      0.5708945
0.55478998      0.54603175      0.52064545     0.56200272      0.63596243
0.88365177      0.96311579
;
var tfplac; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.23485445
 0.23599867      0.23552658      0.23597783 0.23260414      0.22617245
0.21087118      0.19861263      0.19518601     0.2093886       0.24472491
0.24102078      0.26885901
;
var tfpmen; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.1773513
 0.17863369      0.17902194      0.17994256 0.17599616      0.1688019
0.15269466      0.13800441      0.13742846     0.15013 0.16938041
0.1666316       0.1884702
;
var tfprus; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.25550934
 0.25598589      0.25468463      0.25535451 0.25085884      0.24251532
0.23429214      0.22333872      0.23918689     0.25441269      0.2796369
0.43555429      0.24294775
;
var tfpssa; periods 1 2 3 4 5 6 7 8 9 10 11 12 13:33; values    0.072830787
 0.07363556      0.074326737     0.075369592 0.074512119     0.072238779
0.065726662     0.05982857      0.058836395     0.057469768     0.064671243
 0.057775313     0.05747975
;



end;


options_.maxit = 200;
options_.slowc = 1;



simul(periods=80);



%dynatype(bslur3dyna);

rplot tauadv;
/*rplot zu2adv;
rplot labadv;
rplot bsadv;
*/