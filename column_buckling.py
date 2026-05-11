import numpy as np
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa
    
    Return: העומס P בניוטון (float)
    """
    # חישוב עומס קריסה של אוילר (לפי הקדמה הנדסית) כחסם עליון
    P_euler = (np.pi**2 * E * A * r**2) / (L**2)
    
    # הגדרת פונקציית ההפרש לאיפוס
    def f(P):
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        sec_theta = 1.0 / np.cos(theta)
        sigma_max = (P / A) * (1 + (e * c / r**2) * sec_theta)
        return sigma_max - sigma_allow
        
    # מציאת העומס הקריטי באמצעות bisection מ-0 ועד כמעט עומס אוילר
    P_critical = bisect(f, 1e-6, 0.9999 * P_euler)
    
    return float(P_critical)
