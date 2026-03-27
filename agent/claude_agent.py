"""Metabolic engineering agent built on the Claude Agent SDK."""

import sys
from pathlib import Path

import anyio
from dotenv import load_dotenv

# Load ANTHROPIC_API_KEY from .env
load_dotenv(Path(__file__).resolve().parent.parent / ".env")

# Ensure project root is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from claude_agent_sdk import (
    AssistantMessage,
    ClaudeAgentOptions,
    ClaudeSDKClient,
    ResultMessage,
    TextBlock,
    create_sdk_mcp_server,
    tool,
)
from tools.fba_tool import (
    simulate_knockout,
    change_media_and_simulate,
    add_heterologous_reaction,
    simulate_overexpression,
)
from tools.dna_tool import optimize_sequence

MODEL_NAME_DESC = (
    "Name of the metabolic model to use. "
    "Available: 'e_coli_core' (E. coli), 'iMM904' (S. cerevisiae)."
)

# ── FBA tools ────────────────────────────────────────────────────────────────

@tool(
    "simulate_knockout",
    "Run FBA after knocking out one or more genes. "
    "Returns predicted biomass growth rate and flux through a target reaction.",
    {
        "model_name": str,
        "target_reaction": str,
        "knockouts": list[str],
    },
)
async def knockout_tool(args: dict):
    result = simulate_knockout(
        model_name=args["model_name"],
        target_reaction=args["target_reaction"],
        knockouts=args["knockouts"],
    )
    return {"content": [{"type": "text", "text": result}]}


@tool(
    "change_media_and_simulate",
    "Change the growth medium (carbon source and oxygen availability) "
    "then run FBA. Use oxygen_lower_bound=0 for anaerobic, "
    "-1000 for fully aerobic.",
    {
        "model_name": str,
        "carbon_source": str,
        "oxygen_lower_bound": float,
        "target_reaction": str,
    },
)
async def media_tool(args: dict):
    result = change_media_and_simulate(
        model_name=args["model_name"],
        carbon_source=args["carbon_source"],
        oxygen_lower_bound=args["oxygen_lower_bound"],
        target_reaction=args["target_reaction"],
    )
    return {"content": [{"type": "text", "text": result}]}


@tool(
    "add_heterologous_reaction",
    "Add a heterologous (foreign) reaction to the metabolic model "
    "and compute its maximum theoretical flux. "
    "reaction_string uses COBRApy format, e.g. 'pyr_c --> ac_c'.",
    {
        "model_name": str,
        "reaction_id": str,
        "reaction_string": str,
    },
)
async def heterologous_tool(args: dict):
    result = add_heterologous_reaction(
        model_name=args["model_name"],
        reaction_id=args["reaction_id"],
        reaction_string=args["reaction_string"],
    )
    return {"content": [{"type": "text", "text": result}]}


@tool(
    "simulate_overexpression",
    "Simulate gene overexpression by forcing a minimum flux on a reaction, "
    "then observe how biomass growth is affected. "
    "Returns the biomass change to assess metabolic burden.",
    {
        "model_name": str,
        "target_reaction": str,
        "forced_lower_bound": float,
    },
)
async def overexpression_tool(args: dict):
    result = simulate_overexpression(
        model_name=args["model_name"],
        target_reaction=args["target_reaction"],
        forced_lower_bound=args["forced_lower_bound"],
    )
    return {"content": [{"type": "text", "text": result}]}


# ── DNA tool ─────────────────────────────────────────────────────────────────

@tool(
    "optimize_sequence",
    "Codon-optimize a DNA coding sequence for the target host organism. "
    "Preserves the encoded protein, maximises CAI, and removes BsaI sites "
    "to make the sequence Golden-Gate compatible.",
    {
        "sequence": str,
    },
)
async def dna_tool(args: dict):
    result = optimize_sequence(sequence=args["sequence"])
    return {"content": [{"type": "text", "text": result}]}


# ── MCP server wrapping all tools ───────────────────────────────────────────
mcp_server = create_sdk_mcp_server(
    "synbio-tools",
    tools=[
        knockout_tool,
        media_tool,
        heterologous_tool,
        overexpression_tool,
        dna_tool,
    ],
)

SYSTEM_PROMPT = (
    "你是一个顶级的代谢工程专家。"
    "你现在拥有了修改培养基、导入外源反应、模拟过表达和基因敲除的全套能力。"
    "请根据用户需求组合使用这些工具。"
    "可用的代谢模型：'e_coli_core'（大肠杆菌核心模型）和 'iMM904'（酿酒酵母基因组规模模型）。"
    "请根据用户讨论的物种自动选择对应模型。"
    "拿到数据后，再向用户解释机理并给出工程建议。"
)


# ── Core agent function ─────────────────────────────────────────────────────
async def run_agent(user_input: str) -> str:
    """Send a prompt to the agent and return the final text response."""
    options = ClaudeAgentOptions(
        model="claude-sonnet-4-5",
        system_prompt=SYSTEM_PROMPT,
        mcp_servers={"synbio": mcp_server},
        max_turns=20,
    )

    result_text = ""
    async with ClaudeSDKClient(options=options) as client:
        await client.query(user_input)
        async for message in client.receive_response():
            if isinstance(message, ResultMessage):
                result_text = message.result
            elif isinstance(message, AssistantMessage):
                for block in message.content:
                    if isinstance(block, TextBlock):
                        print(f"[Agent] {block.text}")

    return result_text


# ── Interactive CLI ─────────────────────────────────────────────────────────
async def main():
    print("Metabolic Engineering Agent (Claude Agent SDK)")
    print("Available models: e_coli_core, iMM904")
    print("Type 'quit' to exit.")
    print("-" * 50)
    while True:
        try:
            user_msg = input("\nYou: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nBye!")
            break
        if not user_msg or user_msg.lower() in ("quit", "exit"):
            print("Bye!")
            break
        answer = await run_agent(user_msg)
        print(f"\nAgent:\n{answer}")


if __name__ == "__main__":
    anyio.run(main)
