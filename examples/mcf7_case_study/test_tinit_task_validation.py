"""
MCF7 tINIT Task Validation Test Script

Human-GEM + CCLE MCF7 expression data를 사용하여 tINIT reconstruction 완료 후
필수 배지(medium) 조건에 포함된 반응들의 metabolic task 수행 여부를 즉시 확인한다.

Usage:
    python test_tinit_task_validation.py
    python test_tinit_task_validation.py --sample ACH-000019
    python test_tinit_task_validation.py --no-verbose
"""
import os
import argparse

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH  = os.path.join(BASE_DIR, "data", "Human-GEM.xml")
EXPR_PATH   = os.path.join(BASE_DIR, "data", "CCLE_expression_ACH-000019_scores.csv")
TASK_PATH   = os.path.join(BASE_DIR, "data", "metabolicTasks_Essential.xlsx")
RESULTS_DIR = os.path.join(BASE_DIR, "results", "task_validation")


def main(args):
    from troppo import TINITReconstructor

    # --- 사전 파일 확인 ---
    for label, path in [("Model", MODEL_PATH), ("Expression", EXPR_PATH), ("Tasks", TASK_PATH)]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} file not found: {path}")

    os.makedirs(RESULTS_DIR, exist_ok=True)

    # --- 1. Reconstructor 초기화 ---
    recon = TINITReconstructor(
        model_path=MODEL_PATH,
        expression_path=EXPR_PATH,
        solver=args.solver,
        verbose=args.verbose,
    )

    # --- 2. Reconstruction + Task Validation 실행 ---
    samples = [args.sample] if args.sample else None
    results = recon.run(
        samples=samples,
        gap_fill=False,         # 빠른 테스트를 위해 gap-fill 생략
        validate_tasks=True,
        task_file=TASK_PATH,
    )

    # --- 3. 결과 저장 및 출력 ---
    task_matrix = results.task_matrix
    if task_matrix is not None and not task_matrix.empty:
        task_csv = os.path.join(RESULTS_DIR, "task_validation_results.csv")
        task_matrix.to_csv(task_csv)
        print(f"\nTask matrix saved  : {task_csv}")

        summary = results.task_summary()
        summary_csv = os.path.join(RESULTS_DIR, "task_validation_summary.csv")
        summary.to_csv(summary_csv)
        print(f"Summary saved      : {summary_csv}")
        print("\n=== Task Validation Summary ===")
        print(summary.to_string())
    else:
        print("No task validation results available.")

    # --- 4. Medium-related tasks 출력 ---
    med_tasks = results.medium_related_tasks
    if med_tasks:
        print(f"\nMedium-related tasks ({len(med_tasks)}):")
        for t in med_tasks:
            print(f"  - {t}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="tINIT reconstruction with immediate metabolic task validation (MCF7)"
    )
    parser.add_argument(
        "--sample", type=str, default=None,
        help="Single sample name to test (default: first sample in expression data)"
    )
    parser.add_argument(
        "--solver", type=str, default="CPLEX",
        choices=["CPLEX", "GUROBI"],
        help="LP/MILP solver (default: CPLEX)"
    )
    parser.add_argument(
        "--no-verbose", dest="verbose", action="store_false", default=True,
        help="Suppress verbose output"
    )
    args = parser.parse_args()
    main(args)
