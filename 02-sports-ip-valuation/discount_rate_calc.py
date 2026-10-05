"""
Zheng Qinwen Commercial IP Valuation — Discount Rate Derivation
Author: Peixuan Cai (FINC5001)
"""

def calculate_discount_rate():
    # 1. Revenue-Weighted Risk-Free Rate
    w_global, rf_us = 0.40, 0.0400
    w_domestic, rf_cn, spread_inf = 0.60, 0.0175, 0.0292
    rf_final = (w_global * rf_us) + (w_domestic * (rf_cn + spread_inf))

    # 2. Equity Risk Premium
    erp_mature, lambda_sports = 0.0550, 0.0030
    erp_final = erp_mature + lambda_sports

    # 3. Pure-Play Beta
    beta_unlevered, alpha_specificity = 1.30, 1.423
    beta_final = beta_unlevered * alpha_specificity

    # 4. Specific Risk Premiums
    crp = 0.0080  # Country Risk Premium
    arp = 0.0250  # Athlete Risk Premium (1.5% injury + 1.0% career)

    # 5. Consolidated Discount Rate
    discount_rate = rf_final + (beta_final * erp_final) + crp + arp

    print(f"Revenue-Weighted Rf : {rf_final:.4f} -> {rf_final*100:.2f}%")
    print(f"Adjusted ERP        : {erp_final:.4f} -> {erp_final*100:.2f}%")
    print(f"Pure-Play Beta      : {beta_final:.4f} -> {beta_final:.2f}")
    print(f"Final Discount Rate : {discount_rate:.4f} -> {discount_rate*100:.2f}% (Rounded: {round(discount_rate*100, 1)}%)")

if __name__ == "__main__":
    calculate_discount_rate()
