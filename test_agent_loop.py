"""Mock-based test for the agentic tool-use loop.

Simulates Anthropic API responses to verify:
  1. tool_use → dispatch → tool_result → re-call flow
  2. Multi-tool chaining (FBA then DNA)
  3. end_turn exits the loop and returns text
  4. Real tool functions are executed (cobra / dnachisel)
"""

import copy
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))


# ── Helpers to build fake Anthropic response objects ─────────────────────────
def _text_block(text: str):
    return SimpleNamespace(type="text", text=text)


def _tool_use_block(tool_id: str, name: str, tool_input: dict):
    return SimpleNamespace(type="tool_use", id=tool_id, name=name, input=tool_input)


def _make_response(content: list, stop_reason: str):
    return SimpleNamespace(content=content, stop_reason=stop_reason)


# ── The test ─────────────────────────────────────────────────────────────────
def test_full_agent_loop():
    """Simulate a 3-round conversation:
       Round 1: Claude calls simulate_knockout  (stop_reason=tool_use)
       Round 2: Claude calls optimize_sequence   (stop_reason=tool_use)
       Round 3: Claude returns final text         (stop_reason=end_turn)
    """

    responses = [
        # Round 1 – FBA tool call
        _make_response(
            content=[
                _text_block("Let me first check the metabolic flux."),
                _tool_use_block(
                    "call_001",
                    "simulate_knockout",
                    {"target_reaction": "PFK", "knockouts": ["b0721"]},
                ),
            ],
            stop_reason="tool_use",
        ),
        # Round 2 – DNA tool call
        _make_response(
            content=[
                _text_block("Now let me optimize the gene sequence."),
                _tool_use_block(
                    "call_002",
                    "optimize_sequence",
                    {"sequence": "ATGAAAGCAATTTTCGTACTGAAAGATTTCGCCAATAAATAA"},
                ),
            ],
            stop_reason="tool_use",
        ),
        # Round 3 – Final answer
        _make_response(
            content=[_text_block("Based on FBA and codon optimisation results, here is my analysis.")],
            stop_reason="end_turn",
        ),
    ]

    # Capture a deep-copy of the messages kwarg at each call
    captured_messages = []
    call_count = 0

    def fake_create(**kwargs):
        nonlocal call_count
        captured_messages.append(copy.deepcopy(kwargs["messages"]))
        resp = responses[call_count]
        call_count += 1
        return resp

    with patch("agent.claude_agent.client") as mock_client:
        mock_client.messages.create.side_effect = fake_create

        from agent.claude_agent import run_agent
        result = run_agent("敲除 b0721 后 PFK 通量如何？请优化相关基因序列。")

    # ── Assertions ───────────────────────────────────────────────────────────
    assert call_count == 3, f"Expected 3 API calls, got {call_count}"

    # Round 1: just the user message
    msgs_r1 = captured_messages[0]
    assert len(msgs_r1) == 1
    assert msgs_r1[0]["role"] == "user"
    print("[Round 1] 1 message (user) → tool_use(simulate_knockout)")

    # Round 2: user + assistant(tool_use) + user(tool_result)
    msgs_r2 = captured_messages[1]
    assert len(msgs_r2) == 3
    assert msgs_r2[1]["role"] == "assistant"
    assert msgs_r2[2]["role"] == "user"
    fba_result = msgs_r2[2]["content"][0]["content"]
    assert "Biomass flux" in fba_result, f"FBA result missing biomass: {fba_result}"
    assert "PFK" in fba_result
    print(f"[Round 2] 3 messages → tool_use(optimize_sequence)")
    print(f"  FBA tool_result: {fba_result}")

    # Round 3: 5 messages total
    msgs_r3 = captured_messages[2]
    assert len(msgs_r3) == 5
    dna_result = msgs_r3[4]["content"][0]["content"]
    assert "Optimized sequence" in dna_result
    print(f"[Round 3] 5 messages → end_turn")
    print(f"  DNA tool_result: {dna_result}")

    # Final returned text
    assert "analysis" in result.lower()
    print(f"\n[Final output] {result}")

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    test_full_agent_loop()
