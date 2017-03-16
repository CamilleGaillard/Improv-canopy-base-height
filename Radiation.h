// ----------------------------------------------------------------------------

// Berechne den relativen inversen Sonne-Erde-Abstand in Abhängigkeit 
// vom Jahrestag (ca. 0,97<dr<1.03)
double FOAdr(int DayofYear);

// Berechne den Deklinationswinkel der Sonne in Abh. vom Jahrestag
// Deklination : Winkel zwischen Objekt und Aquatorebene
// (ca. -0,4 < delta < 0,4, entspricht ca. -22 bis 22 Grad)
double FOAdelta(int DayofYear);

// Berechne den Winkel zur Sonne bei Sonnenuntergang in Abh. von
// der geogr. Breite psi und der Deklination delta
double FOAOmega(double psi, double delta);

// Berechnet die tageweise extraterrestrische Sonneneinstrahlung in Abh. von
// der Sonnenkonstante Gsc=0.082, der inversen relativen Sonnendistanz dr,
// dem Sonnenuntergangswinkel omega, der geogr. Breite sowie der Deklination delta.
double FOARa(double dr, double omega, double psi, double delta);

// Berechnet die Sonnenstunden pro Tag in Abhängigkeit vom Sonnenuntergangswinkel
double FOADaylighthours(double omega);

// Berechnet eine Schätzung für die Sonnenstrahlung, dabei ist Ra die extraterr.
// Strahlung, Angstrong* sind Regressionsparameter, n ist die aktuelle Sonnenschein-
// dauer am Tag und dlh ist die maximal mögliche Sonnenscheindauer am Tag.
double FOARs(double Ra, double n, double dlh);

// Berechnet eine Schätzung für die Sonnenstrahlung bei klarem Himmel, d.h., es 
// gibt keine Gewichtung mittls Sonnenscheindauer am Tag (s. Rs).
double FOARso(double Ra);

double FOARsoII(double Ra, double z);

// Berechnetdie Netto-Sonneneinstrahlung abhängig von der eingehenden Sonnen-
// einstrahlung und dem ALBEDO.
double FOARns(double Rs);

// Berechnet die langwellige Strahlung abhängig von der Sonnenstrahlung Rs, der 
// maximal möglichen Strahlung Rso, von Max./Min.-Temperatur, von dem Dampfdruck eA
// und der Luftdichte sigma. 
// double FOARnl(double Rs, double Rso, double Tmax, double Tmin, double eA, double sigma);
double FOARnl(double Rs, double Rso, double Tmax, double Tmin, double eA);

// Berechnet die Netto-Strahlung in Abhängigkeit von diversen Daten: Ort, Temperatur,
// Dampfdruck, ALBEDO, Zeit. In der Funktion werden die übrigen definierten Funktionen
// aufgerufen, um wichtige Zwischenwerte zu berechnen.
void   GetNetRadiation(double latitude, double psunshine, 
					   double Tmax, double Tmin, double eA, 
					   int day, double *Rn, double *SunHrs);


double GetPARRadiation( double latitude, double psunshine, int day );
					  
