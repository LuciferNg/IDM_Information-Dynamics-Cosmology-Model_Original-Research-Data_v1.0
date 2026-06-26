# Hermes alert: Multi-test batch for IDM cosmology

Loading: userskills-orchestrator-workflow, userskills-response-conventions

This is a large computational task with 4 parallel sub-tests (BAO, growth rate, likelihood contours, Alcock-Paczynski). Use orchestrator-agent workflow: one coordinator + 3 subagents for the independent tests.

Tests are independent and can run in parallel.
