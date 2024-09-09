import torch
import torch.nn.functional as F


class NNRG(torch.nn.Module):
      def __init__(self,i_dim=640, h_dim=16, o_dim=1):
          super(NNRG, self).__init__()
          self.linear1 = torch.nn.Linear(i_dim, h_dim)
          self.linear2 = torch.nn.Linear(h_dim, o_dim)

      def forward(self, x):
          x = F.relu(self.linear1(x))
          y_pred = self.linear2(x)
          return y_pred







